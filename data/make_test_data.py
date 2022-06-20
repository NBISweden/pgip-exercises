#!/usr/bin/env python3
import argparse
import itertools
import multiprocessing
import os
import pathlib
import copy
import re
import json
import subprocess
import tempfile
import numpy as np
from dataclasses import dataclass, field
from random import choices

import cyvcf2
import tskit


def _get_metadata(tablerow):
    md_string = tablerow.metadata.decode()
    if md_string:
        md = json.loads(tablerow.metadata.decode())
    else:
        md = {}
    return md

def _kw_options(kw):
    args = []
    for k, v in kw.items():
        if len(k) == 1:
            args += [f"-{k}"]
        else:
            args += [f"--{k}"]
        args += [str(v)]
    return args


def samtools_faidx(fasta):
    subprocess.run(["samtools", "faidx", fasta, "-o", f"{str(fasta)}.fai"], check=True)

def bwa_index(fasta):
    subprocess.run(["bwa", "index", fasta], check=True)

def gatk_create_sequence_dictionary(fasta):
    d = str(fasta).replace(".fasta", ".dict")
    if os.path.exists(d):
        os.unlink(d)
    subprocess.run(["gatk", "CreateSequenceDictionary", "-R", fasta], check=True)


@dataclass
class SeqRecord:
    id: str
    name: str
    seq: str
    haplotype: int = -1
    node: int = -1
    population: str = None

    @property
    def seqid(self):
        return f"{self.population}-{self.id}-{self.haplotype}"

    @property
    def description(self):
        return (
            f"id:{self.id}, node:{self.node}, name:{self.name}, "
            f"haplotype:{self.haplotype}, "
            f"population:{self.population}"
        )

    def as_array(self):
        return list(self.seq)

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        s = f">{self.id} {self.description}\n"
        n = 72
        chunks = [self.seq[i : i + n] for i in range(0, len(self), n)]
        s += "\n".join(chunks)
        s += "\n"
        return s

    def __repr__(self):
        return str(self)


@dataclass
class Individual:
    uniqueid: int
    id: int
    population: int
    population_name: str = None
    metadata: dict = None
    haplotype: list[SeqRecord] = field(default_factory = list)
    is_reference: bool = False
    readfiles: tuple = ()
    bamfile: str = None
    gvcf: str = None

    def __str__(self):
        return f"{self.population_name}-{self.id}"

    def __repr__(self):
        return (f"Individual(uniqueid={self.uniqueid}, id={self.id}, "
                f"population={self.population}, population_name={self.population_name}, "
                f"metadata={self.metadata})")

    @property
    def tskit_id(self):
        return f"tsk_{self.uniqueid}"

    @property
    def reference(self):
        try:
            return self.haplotype[0]
        except Excetption as e:
            print(e)
            pass

    def make_sequence(self, vcf, reference):
        if self.is_reference:
            return
        fh = cyvcf2.VCF(vcf, samples=self.tskit_id)
        for i, haplotype in enumerate(reference.haplotype):
            seqid = f"{self.tskit_id}-{self.population_name}-{self.id}-{i}"
            sequence = copy.deepcopy(haplotype.as_array())
            for gt in fh:
                if gt.genotypes[0][i] == 1:
                    sequence[gt.POS - 1] = gt.ALT[0]
            record = SeqRecord(
                id=seqid,
                name=seqid,
                seq="".join(sequence),
                haplotype=i,
                population=self.population_name
            )
            self.haplotype.append(record)

    def write_haplotypes(self, path):
        with open(path, "w") as fh:
            for hap in self.haplotype:
                fh.write(str(hap))

    def simulate_reads(self, path, *, output, coverage=10.0, readlength=125, **kw):
        nreads = int((len(self.haplotype[0]) * coverage) / (readlength * 2))
        tmp = tempfile.mkstemp(suffix=".fasta")[1]
        self.write_haplotypes(tmp)
        args = ["iss", "generate", "--output", output,  "--coverage", "uniform",
                "--genomes", tmp, "-z", "--model", "HiSeq", "-n", str(nreads)]
        args += _kw_options(kw)
        subprocess.run(args, check=True)
        os.unlink(tmp)
        os.unlink(f"{output}_coverage.txt")
        self.readfiles = (f"{output}_R1.fastq.gz", f"{output}_R2.fastq.gz")

    def align_reads(self, idxbase, **kw):
        self.bamfile = self.readfiles[0].replace("_R1.fastq.gz", ".bam")
        args = ["bwa", "mem", "-R", f"\'@RG\\tID:{self.tskit_id}\\tSM:{str(self)}\'"]
        args += _kw_options(kw) + [idxbase,  self.readfiles[0], self.readfiles[1]]
        args += ["|", "samtools", "fixmate", "-m", "-", "/dev/stdout"]
        args += ["|", "samtools", "sort", "-"]
        args += ["|", "samtools", "markdup", "-", "/dev/stdout"]
        args += ["|", "samtools", "view", "-h", "-b", "-o", self.bamfile]
        subprocess.run(" ".join(args), check=True, shell=True)
        subprocess.run(["samtools", "index", self.bamfile], check=True)
        for readfile in self.readfiles:
           os.unlink(readfile)
        self.readfiles = ()

    def haplotype_caller(self, idxbase, **kw):
        self.gvcf = self.bamfile.replace(".bam", ".g.vcf")
        args = ["gatk", "HaplotypeCaller", "-I", self.bamfile, "-ERC", "GVCF", "-O", self.gvcf, "-R", idxbase] + _kw_options(kw)
        subprocess.run(args, check=True)


@dataclass
class Population:
    id: int
    name: str
    metadata: dict = field(default_factory = dict)
    individuals: list[Individual] = field(default_factory = list)

    @property
    def size(self):
        return len(self.individuals)



@dataclass
class DemesModel:
    name: str
    path: pathlib.Path
    demesfile: str
    ts: tskit.trees.TreeSequence = None
    populations: list[Population] = field(default_factory=list)
    individuals: list[Individual] = field(default_factory=list)
    vcf_ref: str = None

    def __post_init__(self):
        self._census = 0

    @property
    def census(self):
        return self._census

    @census.setter
    def census(self, census):
        self._census = census

    def add_individuals(self, **kw):
        for pop, count in kw.items():
            if pop in self.popdict.keys():
                population = self.get_population(pop)
            else:
                population = Population(id=len(self.populations), name=pop)
                self.populations.append(population)
            popindex = population.id
            popsize = len(self.populations[popindex].individuals)
            for i in np.arange(popsize, popsize + count):
                ind = Individual(uniqueid=self.census, id=i, population=popindex, population_name=pop)
                self.populations[popindex].individuals.append(ind)
                self.individuals.append(ind)
                self.census += 1

    @property
    def reference(self):
        # loop samples
        for ind in self.individuals:
            if ind.is_reference:
                return ind

    @property
    def index(self):
        """Reference index base name"""
        return self.path / f"{str(self.reference)}.fasta"

    def set_reference(self, individual: int, population: str):
        for pop in self.populations:
            if pop.name == population:
                pop.individuals[individual].is_reference = True

    def make_reference_sequence(self):
        """Create the reference sequence for the reference individual.
        Assume biallelic snps (no back mutations).

        Haplotype 0 is the reference; haplotype 1 the reference with variants.
        """
        if self.ts is None:
            return
        if self.vcf_ref is None:
            return
        dna = ["A", "C", "G", "T"]
        hap0 = choices(dna, k=int(self.ts.sequence_length))
        hap1 = copy.deepcopy(hap0)
        # We must make sure the reference holds the actual reference
        # allele at the variant positions as obtained from the tree
        # sequence
        for gt in cyvcf2.VCF(self.vcf_ref):
            ref = gt.genotypes[0][0]
            hap0[gt.POS - 1] = gt.REF
            hap1[gt.POS - 1] = gt.ALT[0]
        h0 = f"{self.reference.tskit_id}-{self.reference.population_name}-{self.reference.id}-0"
        h1 = f"{self.reference.tskit_id}-{self.reference.population_name}-{self.reference.id}-1"
        self.reference.haplotype = [SeqRecord(id=h0, name=h0, seq="".join(hap0), haplotype=0, population=self.reference.population_name),
                                    SeqRecord(id=h1, name=h1, seq="".join(hap1), haplotype=1, population=self.reference.population_name)]


    def save_reference(self):
        with open(self.index, "w") as fh:
            fh.write(str(self.reference.haplotype[0]))


    def get_population(self, name):
        for pop in self.populations:
            if pop.name == name:
                return pop

    @property
    def popdict(self):
        return dict((p.name, p.size) for p in self.populations)

    # FIXME: should be agnostic to msprime / SLiM
    def simulate(self, *, seqlength, recombination_rate,
                 mutation_rate, anc_seed=42, mut_seed=10):
        tsfile = self._ts_ancestry(seqlength, recombination_rate, anc_seed)
        tsfile = self._ts_mutate(tsfile, mutation_rate, mut_seed)
        self.ts = tskit.load(tsfile)
        os.unlink(tsfile)

    def _ts_ancestry(self, seqlength, recombination_rate, seed):
        populations = [f"{k}:{v}" for k, v in self.popdict.items()]
        tsfile = self.path / f"{self.name}.ts"
        args = ["msp", "ancestry",
                "--demography", self.demesfile,
                "--length", str(seqlength),
                "--recombination-rate", str(recombination_rate),
                "--random-seed", str(seed),
                "-o", tsfile,
            ] + populations

        # Run msp ancestry
        subprocess.run(args, check=True)
        return tsfile

    def _ts_mutate(self, tsfile, mutation_rate, seed):
        tsmutfile = self.path / f"{self.name}.mut.ts"
        args = ["msp", "mutations",
                "-o", tsmutfile,
                "--random-seed", str(seed),
                str(mutation_rate), tsfile]
        subprocess.run(args, check=True)
        os.unlink(tsfile)
        return tsmutfile

    # FIXME: should sync back to Individual?
    def sync_metadata(self):
        """Sync tree sequence and DemesModel metadata"""
        tc = self.ts.dump_tables()
        indrows = list(getattr(tc, "individuals"))
        getattr(tc, "individuals").reset()
        oldname = None
        for i, row in enumerate(indrows):
            ind = self.ts.individuals()[i]
            popindex = self.ts.nodes()[ind.nodes[0]].population
            pop = self.ts.population(popindex)
            popname = pop.metadata["name"]
            if oldname is not None and oldname != popname:
                i = 0
            oldname = popname
            name = f"tsk_{ind.id}_{popname}_{i}"
            md = {
                "id": ind.id,
                "tskit_id": f"tsk_{ind.id}",
                "name": name,
                "popname": popname,
                "popindex": popindex,
                "ind_index": i,
                "is_reference": self.individuals[i].is_reference,
            }
            md["vcfheader"] = (
                f"##<SAMPLE=<ID={md['tskit_id']},Name={md['name']},"
                f"Index={md['id']},Population={md['popname']},"
                f"Description=\"{pop.metadata['description']}\">"
            )
            row = row.replace(metadata=json.dumps(md).encode())
            getattr(tc, "individuals").append(row)
        self.ts = tc.tree_sequence()


    def write_variants(self):
        """Write variants to output vcf.

        If reference vcf is provided flip alleles where necessary.
        WARNING: this means tree sequence states and vcf states
        will differ!

        """
        vcf_tmp_fn = tempfile.mkstemp(suffix=".vcf")[1]
        print(f"Writing ts to vcf {vcf_tmp_fn}")
        with open(vcf_tmp_fn, "w") as fh:
            self.ts.write_vcf(fh)

        vcf_ref = None
        if self.reference is not None:
            vcf_ref = cyvcf2.VCF(vcf_tmp_fn, samples=self.reference.tskit_id)

        vcf_out = self.path / f"{self.name}.vcf.gz"
        print(f"Writing output vcf {vcf_out}")
        samples = [ind.tskit_id for ind in self.individuals]
        vcf_tmp = cyvcf2.VCF(vcf_tmp_fn, samples=samples)
        vcfwriter = cyvcf2.cyvcf2.Writer(str(vcf_out), vcf_tmp, "wz")
        for ind in self.individuals:
            tsk_ind = self.ts.individuals()[ind.uniqueid]
            vcfwriter.add_to_header(_get_metadata(tsk_ind)["vcfheader"])
        if vcf_ref is None:
            vcfwriter.set_samples(vcf_tmp.samples)
            for variant in vcf_tmp:
                vcfwriter.write_record(variant)
        else:
            # in cyvcf2: Pick one individual *haplotype* as reference: this
            # individual should have only 0's, so all calls at a site with a
            # derived allele should be flipped for all individuals.
            for gt, variant in zip(vcf_ref, vcf_tmp):
                ref = gt.genotypes[0][0]
                if ref == 1:
                    # Flip states for all other genotypes
                    for i in range(len(variant.genotypes)):
                        alleles = variant.genotypes[i]
                        variant.genotypes[i] = [1 - alleles[0], 1 - alleles[1], alleles[2]]
                    variant.genotypes = variant.genotypes
                    variant.REF = gt.ALT[0]
                    variant.ALT = [gt.REF]
                vcfwriter.write_record(variant)
        vcfwriter.close()
        subprocess.run(["tabix", "-f", vcf_out], check=True)
        self.vcf_ref = vcf_out

    # FIXME: dump more than ts?
    def dump(self, **kw):
        path = self.path / f"{self.name}.ts"
        print("Dumping tree sequences to", str(path))
        self.ts.dump(path, **kw)

    def write_haplotypes(self):
        for ind in self.individuals:
            path = self.path / f"{str(ind)}"
            ind.write_haplotypes(path)

    def gatk_variant_calling(self, cpus=1, **kw):
        for ind in self.individuals:
            if ind.is_reference:
                continue
            ind.haplotype_caller(str(self.index), **{'native-pair-hmm-threads': cpus})
        gvcf = str(self.path / f"{self.name}.gatk.g.vcf.gz")
        args = ["gatk", "CombineGVCFs", "-O", gvcf, "-R", self.index]
        gvcflist = []
        for ind in self.individuals:
            if not ind.is_reference:
                gvcflist += ["-V", ind.gvcf]
        subprocess.run(args + gvcflist, check=True)
        for ind in self.individuals:
            if ind.is_reference:
                continue
            os.unlink(ind.gvcf)
            idx = f"{ind.gvcf}.idx"
            if os.path.exists(idx):
                os.unlink(idx)
        vcf = str(self.path / f"{self.name}.gatk.vcf.gz")
        args = ["gatk", "GenotypeGVCFs", "-R", self.index, "-V", gvcf, "-O", vcf]
        subprocess.run(args, check=True)
        subprocess.run(["tabix", "-f", vcf], check=True)


def make_ooa(args):
    np.random.seed(52)
    model = DemesModel(name="ooa", path=pathlib.Path("ooa"), demesfile="ooa/ooa_with_outgroup.demes.yaml")
    popdict = {'CHB': 3, 'CEU': 4, 'YRI': 3, 'gorilla': 1,
               'chimpanzee': 1, 'orangutan': 1}
    model.add_individuals(**popdict)
    model.set_reference(population="CEU", individual=3)
    model.simulate(seqlength=5e4, recombination_rate=1e-7, mutation_rate=5e-9)
    model.dump()
    model.sync_metadata()
    model.write_variants()
    model.make_reference_sequence()
    model.save_reference()
    bwa_index(model.index)
    samtools_faidx(model.index)
    gatk_create_sequence_dictionary(model.index)

    for ind in model.individuals:
        if ind.is_reference:
            continue
        ind.make_sequence(model.vcf_ref, model.reference)
        output = model.path / f"{str(ind)}-0"
        coverage = np.clip(np.random.normal(10, 2), 5.0, 15.0)
        ind.simulate_reads(model.path, cpus=args.cpus, output=output, coverage=coverage)
        ind.align_reads(str(model.index), t=args.cpus)

    model.gatk_variant_calling(cpus=args.cpus)



def main():
    parser = argparse.ArgumentParser(description="Make test data")
    parser.add_argument(
        "--force",
        "-F",
        help=("force dataset creation"),
        action="store_true",
        default=False,
    )

    subparsers = parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    ooaparser = subparsers.add_parser("ooa-outgroup", help="Make ooa dataset (neutral)")
    ooaparser.set_defaults(runner=make_ooa)
    ooaparser.add_argument(
        "--cpus", "-p", help="Set number of cpus", default=1, dest="cpus"
    )

    args = parser.parse_args()
    args.runner(args)


if __name__ == "__main__":
    main()
