package com.epam.prorecon;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

public class FileReaderUtils {
    public static String readSequenceFromFastaFile(String fileUrl, String chromosomeName, Integer beginIndex,
                                               Integer endIndex) throws FileNotFoundException {
        final File refFile = new File(fileUrl);
        IndexedFastaSequenceFile sequenceFile = new IndexedFastaSequenceFile(refFile);

        String nucleotideString = sequenceFile.getSubsequenceAt(chromosomeName, beginIndex, endIndex).getBaseString();
        CloserUtil.close(sequenceFile);
        return nucleotideString;
    }

    public static List<VariantContext> readVariantContextsFromVcfFile(String fileUrl, String chromosomeName,
                                                                      Integer beginIndex, Integer endIndex) {
        final File vcfFile = new File(fileUrl);
        VCFFileReader vcfFileReader = new VCFFileReader(vcfFile, false);

        List<VariantContext> variantContexts = vcfFileReader.query(chromosomeName, beginIndex, endIndex).toList();
        CloserUtil.close(vcfFile);
        return variantContexts;
    }
}
