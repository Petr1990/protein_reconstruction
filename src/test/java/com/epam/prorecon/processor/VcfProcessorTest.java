package com.epam.prorecon.processor;

import com.epam.prorecon.FileReaderUtils;
import com.epam.prorecon.entity.ExonPosition;
import com.epam.prorecon.entity.WorkflowResult;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.bed.FullBEDFeature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.vcf.VCFCodec;
import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

public class VcfProcessorTest {
    private static final Logger LOGGER = LoggerFactory.getLogger(VcfProcessorTest.class);

    @Test
    public void nucliotideStringAfterApplyingVcfShouldBeEqualToPrecalculatedConstant() throws IOException, CloneNotSupportedException {
        String fastaFileSubSequence = FileReaderUtils.readSequenceFromFastaFile(
                this.getClass().getClassLoader().getResource("dmel-all-chromosome-r606.fasta"/*"Ref1.fasta"*/)
                        .getPath(), "X", 12584385, 12592193);

        URL vcfFileUrl = this.getClass().getClassLoader().getResource("agnX1.model.2.snp-indels.vcf"/*"mutatsii.vcf"*/);
        Assert.assertNotNull(vcfFileUrl);
        File vcfFile = new File(vcfFileUrl.getPath());

        URL gtfFileUrl = this.getClass().getClassLoader().getResource("dmel-all-r6.06.LIMK1.gtf");
        Assert.assertNotNull(gtfFileUrl);

        File vcfIndexFile = new File(vcfFileUrl.getPath() + ".Idx");
        Index idx = IndexFactory.createIndex(vcfFile, new VCFCodec(), IndexFactory.IndexType.LINEAR);

        LittleEndianOutputStream stream = null;
        stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(vcfIndexFile)));
        idx.write(stream);
        stream.close();

        CloserUtil.close(vcfFile);
        CloserUtil.close(vcfIndexFile);

        List<ExonPosition> exonPositions = FileReaderUtils.readGffFile(gtfFileUrl.getPath());
        VcfProcessor vcfProcessor = new VcfProcessor(FileReaderUtils.readVariantContextsFromVcfFile(
                vcfFileUrl.getPath(), "X", 12584385, 12592193));

        vcfProcessor.addReferenceResult(fastaFileSubSequence, exonPositions, 12584385);
        vcfProcessor.process(fastaFileSubSequence, fastaFileSubSequence, 12584385, 12584385, exonPositions,
                new ArrayList<>(exonPositions));
        PrintWriter out = new PrintWriter("C:\\Users\\user\\Downloads\\2_version_results_1_petr.txt");
        for (WorkflowResult workflowResult : vcfProcessor.getWorkflowResults()) {
//            LOGGER.warn(currString);
            out.println(fastaFileSubSequence);
            out.println(workflowResult.getAppliedMutationsDnaString());
            out.println(workflowResult.getMrnaString());
            out.println(workflowResult.getMrnaWithoutIntronsString());
            out.println(workflowResult.getProteinString());
            out.println();
        }
        out.close();
        LOGGER.warn(fastaFileSubSequence);
    }
}
