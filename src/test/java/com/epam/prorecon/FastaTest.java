package com.epam.prorecon;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

public class FastaTest {
    private final static String FIRST_FASTA_SEQUANCE_BEGIN =
            "TAATTAAAACAGATCCTGAGAAAATTTCCACAATTATGAAGTATCCGATTCCACAAAACATTAGAGA" +
            "GCTTCGAAGTTTTCTAGGCCTCACCGGCTACTACCGTAAATT";

    @Test
    public void fastaFileStartShouldEqualToReferenceString() throws IOException {
        final File fasta = new File(this.getClass().getClassLoader().getResource("dmel-all-chromosome-r606.fasta").getPath());
        final FastaSequenceFile fastaReader = new FastaSequenceFile(fasta, true);
        final ReferenceSequence referenceSequence = fastaReader.nextSequence();
        Assert.assertEquals(referenceSequence.getName(), "211000022279114");
        Assert.assertEquals(StringUtil.bytesToString(referenceSequence.getBases()).substring(0,
                FIRST_FASTA_SEQUANCE_BEGIN.length()), FIRST_FASTA_SEQUANCE_BEGIN);
    }
}
