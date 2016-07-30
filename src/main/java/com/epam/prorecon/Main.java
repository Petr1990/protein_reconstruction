package com.epam.prorecon;

import com.epam.prorecon.entity.ExonPosition;
import com.epam.prorecon.processor.SequenceProcessor;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.vcf.VCFCodec;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class Main {
    public static void main(String[] args) throws IOException, CloneNotSupportedException {
        String fastaFileSubSequence = FileReaderUtils.readSequenceFromFastaFile(
                new File(args[0]).getPath(), args[3], Integer.valueOf(args[4]), Integer.valueOf(args[5]));

        File vcfFile = new File(args[1]);

        File gtfFileUrl = new File(args[2]);

        File vcfIndexFile = new File(vcfFile.getPath() + ".Idx");
        Index idx = IndexFactory.createIndex(vcfFile, new VCFCodec(), IndexFactory.IndexType.LINEAR);

        LittleEndianOutputStream stream = null;
        stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(vcfIndexFile)));
        idx.write(stream);
        stream.close();

        CloserUtil.close(vcfFile);
        CloserUtil.close(vcfIndexFile);

        List<ExonPosition> exonPositions = FileReaderUtils.readGffFile(gtfFileUrl.getPath());
        SequenceProcessor sequenceProcessor = new SequenceProcessor(FileReaderUtils.readVariantContextsFromVcfFile(
                vcfFile.getPath(), args[3], Integer.valueOf(args[4]), Integer.valueOf(args[5])),
                Integer.valueOf(args[4]));

        sequenceProcessor.addReferenceResult(fastaFileSubSequence, exonPositions, Integer.valueOf(args[4]));
        sequenceProcessor.process(fastaFileSubSequence, fastaFileSubSequence, Integer.valueOf(args[4]),
                Integer.valueOf(args[4]), exonPositions, exonPositions.stream().map(ExonPosition::new).
                        collect(Collectors.toList()));
        PrintWriter protreinsFile = new PrintWriter("proteins.fasta");
        PrintWriter totalFile = new PrintWriter("total.txt");
        List<String> proteins = new ArrayList<>();
        totalFile.println(fastaFileSubSequence);
        for (int i = 0; i < sequenceProcessor.getWorkflowResults().size(); i++) {
            totalFile.println(sequenceProcessor.getWorkflowResults().get(i).getAppliedMutationsDnaString());
            totalFile.println(sequenceProcessor.getWorkflowResults().get(i).getProteinString());
            totalFile.println();

            if (!sequenceProcessor.getWorkflowResults().get(i).getProteinString().equals("")
                    && !proteins.contains(sequenceProcessor.getWorkflowResults().get(i).getProteinString())) {
                protreinsFile.println(">" + i);
                proteins.add(sequenceProcessor.getWorkflowResults().get(i).getProteinString());
                protreinsFile.println(sequenceProcessor.getWorkflowResults().get(i).getProteinString());
            }
        }
        protreinsFile.close();
        totalFile.close();
    }
}
