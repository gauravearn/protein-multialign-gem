#!/usr/bin/ruby
=begin
 Author: Gaurav Sablok
 Universitat Potsdam
 Date: 2024-4-5
 a ruby class for analyzing the single and the multigff alignments
 coming from the genome annotations. It has the support for both the
 gff and the paf alignments. Also has the support for building the graphs
 for the graphql.
=end
gem install bio
gem install parse_fasta
require 'bio'
require 'parse_fasta'

class SINGLEMultiGFF

  def extractpafSingle(fastafile, gffalignment)
      fasta = File.read(fastafile).readlines
      alignment = File.read(gffalignment).readlines
      sequenceids = []
      sequencechar = []
      ParseFasta::SeqFile.open(fasta).each_record do | iter |
      sequenceids.push(iter.header.to_s)
      sequencechar.push(iter.seq.to_s)
      end
      gffread = File.open(gffalignment).readlines
      gffreadmain = gffread[2..gffread.length]
      mRNAidlocals = []
      mRNAid_start_positions = []
      mRNAid_stop_position = []
      mRNAgffids = []
      for i in 0..gffreadmain.length
        mRNAidslocal.push(gffreadmain[i].to_s.strip.split[2]) \
                                  if gffread[i].to_s.strip.split[2] == "mRNA"
      end
      for i in 0..gffreadmain.length
        mRNAids_start_position.push(gffreadmain[i].to_s.strip.split[3]) \
                                  if gffread[i].to_s.strip.split[2] == "mRNA"
      end
      for i in 0..gffreadmain.length
        mRNAids_stop_position.push(gffreadmain[i].to_s.strip.split[4]) \
                                 if gffread[i].to_s.strip.split[2] == "mRNA"
      end
      for i in 0..gffreadmain.length
        mRNAgffids.push(gffreadmain[i].to_s.strip.split[0]) \
                                 if gffread[i].to_s.strip.split[2] == "mRNA"
      end
      mRNAsequenceiter = {}
      mRNAsequencesplice = []
      for i in 0..mRNAid_start_positions.length
         mRNsequenceiter[sequenceids[i]] = [mRNAid_start_positions[i], mRNAid_stop_positions]
      end
      for i in 0..mRNAid_start_positions.length
        mRNAsequencesplice[[sequenceids[i],mRNAidslocal[i]] = sequencechar[i].to_s.slice(mRNAid_start_positons[i], mRNAid_end_positions[i])
              if sequenceids[i] == mRNAgffids[i]
      end
  end

  def extractpafMulti(fastafile, multigff, outputfasta)
      fasta = File.read(fastafile).readlines
      alignment = File.read(gffalignment).readlines
      sequenceids = []
      sequencechar = []
      ParseFasta::SeqFile.open(fasta).each_record do | iter |
      sequenceids.push(iter.header.to_s)
      sequencechar.push(iter.seq.to_s)
      end
      multigffread = File.open(multigff, "r").readlines
      multigffmRNA = []
      multigffmRNAstart = []
      multigffmRNAend = []
      multigffseqids = []
      for i in 0..multigffread.length
         multigffmRNA.push(multigffread[i].to_s.strip.split[2]) \
                       if multigffread[i].to_s.strip.split[2] == "mRNA"
      end
      for i in 0..multigffread.length
        multigffmRNAstart.push(multigffread[i].to_s.strip.split[3]) \
                           if multigffread[i].to_s.strip.split[2] == "mRNA"
      end
      for i in 0..multigffread.length
        multigffmRNAend.push(multigffread[i].to_s.strip.split[4]) \
                                if multigffread[i].to_s.strip.split[2] == "mRNA"
      end
      for i in 0..multigffread.length
        multigffseqids.push(multigffread[i].to_s.strip.split[0]) \
                          if multigffread[i].to_s.strip.split[2] == "mRNA"
      end
      multifastaiterinfo = {}
      multifastaseqsplice = {}
      multifastaupstream = []
      for i in 0..multigffseqids.length
         multifastaiterinfo[multigffseqids[i]] = [multigffmRNAstart[i], multigffmRNAend[i]]
      end
      for i in 0..multigffmRNAstart.length
        multifastseqsplice[sequenceids[i]] = sequencechar[i].slice(multigffmRNAstart[i].to_i, multigffmRNAend[i].to_i) \
                       if sequenceids[i].to_s == multigffseqids[i].to_s
      end
  end
  def extractUpstreamDownstream(fastafile, multigff, upstream, downstream)
      fasta = File.read(fastafile).readlines
      alignment = File.read(gffalignment).readlines
      sequenceids = []
      sequencechar = []
      ParseFasta::SeqFile.open(fasta).each_record do | iter |
      sequenceids.push(iter.header.to_s)
      sequencechar.push(iter.seq.to_s)
      end
      multigffread = File.open(multigff, "r").readlines
      multigffmRNA = []
      multigffmRNAstart = []
      multigffmRNAend = []
      multigffseqids = []
      for i in 0..multigffread.length
         multigffmRNA.push(multigffread[i].to_s.strip.split[2]) \
                       if multigffread[i].to_s.strip.split[2] == "mRNA"
      end
      for i in 0..multigffread.length
        multigffmRNAstart.push(multigffread[i].to_s.strip.split[3]) \
                           if multigffread[i].to_s.strip.split[2] == "mRNA"
      end
      for i in 0..multigffread.length
        multigffmRNAend.push(multigffread[i].to_s.strip.split[4]) \
                                if multigffread[i].to_s.strip.split[2] == "mRNA"
      end
      for i in 0..multigffread.length
        multigffseqids.push(multigffread[i].to_s.strip.split[0]) \
                          if multigffread[i].to_s.strip.split[2] == "mRNA"
      end
      multifastaiterinfo = {}
      multifastaseqsplice = {}
      multifastaupstream = {}
      multifastadownstream = {}
      for i in 0..multigffseqids.length
          multifastaiterinfo[multigffseqids[i]] = [multigffmRNAstart[i], multigffmRNAend[i]]
      end
      for i in 0..multigffmRNAstart.length
          multifastaseqsplice[sequenceids[i]] = sequencechar[i].slice(multigffmRNAstart[i].to_i, multigffmRNAend[i].to_i) \
                       if sequenceids[i].to_s == multigffseqids[i].to_s
      end
      for i in 0..multigffmRNAstart.length
          multifastaupstream[sequenceids[i]] = sequencechar[i].slice(multigffmRNAstart[i].to_i+upstream.to_i, multigffmRNAend[i].to_i) \
                       if sequenceids[i].to_s == multigffseqids[i].to_s
      end
      for i in 0..multigffmRNAstart.length
          multifastadownstream[sequenceids[i]] = sequencechar[i].slice(multiffmRNAstart[i].to_i, multigffmRNAend[i].to_i+downstream.to_i) \
                       if sequenceids[i].to_s == multigffseqids[i].to_s
      end
  end
end
