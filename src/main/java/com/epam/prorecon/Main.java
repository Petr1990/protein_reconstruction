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
                new File("dmel-all-chromosome-r606.fasta").getPath(), "X", 12584385, 12592193);

        File vcfFile = new File("agnX1.model.2.snp-indels.vcf");

        File gtfFileUrl = new File("dmel-all-r6.06.LIMK1.gtf");

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
                vcfFile.getPath(), "X", 12584385, 12592193), 12584385);

        sequenceProcessor.addReferenceResult(fastaFileSubSequence, exonPositions, 12584385);
        sequenceProcessor.process(fastaFileSubSequence, fastaFileSubSequence, 12584385, 12584385, exonPositions,
                exonPositions.stream().map(ExonPosition::new).collect(Collectors.toList()));
        PrintWriter protreinsFile = new PrintWriter("proteins.fasta");
        PrintWriter totalFile = new PrintWriter("total.txt");
        List<String> proteins = new ArrayList<>();
        totalFile.println(fastaFileSubSequence);
        for (int i = 0; i < sequenceProcessor.getWorkflowResults().size(); i++) {
            totalFile.println(sequenceProcessor.getWorkflowResults().get(i).getAppliedMutationsDnaString());
//            totalFile.println(sequenceProcessor.getWorkflowResults().get(i).getMrnaString());
//            totalFile.println(sequenceProcessor.getWorkflowResults().get(i).getMrnaWithoutIntronsString());
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
