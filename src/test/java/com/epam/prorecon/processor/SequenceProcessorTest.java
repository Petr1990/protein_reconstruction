package com.epam.prorecon.processor;

import com.epam.prorecon.FileReaderUtils;
import com.epam.prorecon.entity.ExonPosition;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.vcf.VCFCodec;
import org.junit.Assert;
import org.junit.Test;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class SequenceProcessorTest {
    @Test
    public void nucliotideStringAfterApplyingVcfShouldBeEqualToPrecalculatedConstant() throws IOException, CloneNotSupportedException {
        String fastaFileSubSequence = FileReaderUtils.readSequenceFromFastaFile(
                this.getClass().getClassLoader().getResource("dmel-all-chromosome-r606.fasta")
                        .getPath(), "X", 12584385, 12592193);

        URL vcfFileUrl = this.getClass().getClassLoader().getResource("agnX1.model.2.snp-indels.vcf");
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
        SequenceProcessor sequenceProcessor = new SequenceProcessor(FileReaderUtils.readVariantContextsFromVcfFile(
                vcfFileUrl.getPath(), "X", 12584385, 12592193), 12584385);

        sequenceProcessor.addReferenceResult(fastaFileSubSequence, exonPositions, 12584385);
        sequenceProcessor.process(fastaFileSubSequence, fastaFileSubSequence, 12584385, 12584385, exonPositions,
                exonPositions.stream().map(ExonPosition::new).collect(Collectors.toList()));
        List<String> proteins = new ArrayList<>();
        for (int i = 0; i < sequenceProcessor.getWorkflowResults().size(); i++) {
            if (!sequenceProcessor.getWorkflowResults().get(i).getProteinString().equals("")
                    && !proteins.contains(sequenceProcessor.getWorkflowResults().get(i).getProteinString())) {
                proteins.add(sequenceProcessor.getWorkflowResults().get(i).getProteinString());
            }
        }

        Assert.assertEquals(32769, sequenceProcessor.getWorkflowResults().size());
        Assert.assertEquals(19, proteins.size());
    }
}
