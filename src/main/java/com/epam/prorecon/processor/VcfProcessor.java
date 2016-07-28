package com.epam.prorecon.processor;

import com.epam.prorecon.entity.ExonPosition;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class VcfProcessor extends AbstractProcessor {

    private List<VariantContext> leftVcfList;
    private List<ExonPosition> exonPositions;
    private static Set<String> possibleFinalStrings = new HashSet<String>();

    public VcfProcessor(List<VariantContext> leftVcfList, List<ExonPosition> exonPositions) {
        this.leftVcfList = new ArrayList<>(leftVcfList);
        this.exonPositions = new ArrayList<>(exonPositions);
    }

    public void process(String motherNucleotideString, String fatherNucleotideString, int motherStringShift,
                          int fatherStringShift)
            throws CloneNotSupportedException {
        StringBuilder motherNucleotideStringBuffer = new StringBuilder(motherNucleotideString);
        StringBuilder fatherNucleotideStringBuffer = new StringBuilder(fatherNucleotideString);

        while (leftVcfList.size() > 0) {
            VariantContext variantContext = leftVcfList.remove(0);

            if (variantContext.getGenotypes().get(0).isHet()) {
                if (variantContext.getGenotypes().get(0).isPhased()) {
                    String genotypeString = variantContext.getGenotypes().get(0).getGenotypeString(false);
                    if (genotypeString.indexOf("*") < genotypeString.indexOf("|")) {
                        fatherStringShift = applyMutation(variantContext, fatherNucleotideStringBuffer,
                                fatherStringShift);
                    } else {
                        motherStringShift = applyMutation(variantContext, motherNucleotideStringBuffer,
                                motherStringShift);
                    }
                } else {
                    StringBuilder fatherNucleotideStringBufferForCopy = new StringBuilder(fatherNucleotideStringBuffer);
                    int fatherStringShiftForCopy = applyMutation(variantContext, fatherNucleotideStringBufferForCopy,
                            fatherStringShift);

                    VcfProcessor vcfProcessorWithoutCurrentMutation = new VcfProcessor(leftVcfList, exonPositions);
                    vcfProcessorWithoutCurrentMutation.process(motherNucleotideStringBuffer.toString(),
                            fatherNucleotideStringBufferForCopy.toString(), motherStringShift,
                            fatherStringShiftForCopy);

                    motherStringShift = applyMutation(variantContext, motherNucleotideStringBuffer, motherStringShift);
                }
            } else {
                motherStringShift = applyMutation(variantContext, motherNucleotideStringBuffer, motherStringShift);
                fatherStringShift = applyMutation(variantContext, fatherNucleotideStringBuffer, fatherStringShift);
            }
        }

        possibleFinalStrings.add(motherNucleotideStringBuffer.toString());
        possibleFinalStrings.add(fatherNucleotideStringBuffer.toString());
    }

    private static int applyMutation(VariantContext variantContext, StringBuilder nucleotideStringBuilder, int shift) {
        nucleotideStringBuilder.replace(variantContext.getStart() - shift, variantContext.getEnd() - shift + 1,
                variantContext.getAlleles().get(1).getBaseString());
        return shift + variantContext.getAlleles().get(0).getBaseString().length()
                - variantContext.getAlleles().get(1).getBaseString().length();
    }

    public Set<String> getPossibleFinalStrings() {
        return possibleFinalStrings;
    }
}
