package com.epam.prorecon.processor;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class VcfProcessor extends AbstractProcessor {

    private List<VariantContext> leftVcfList;
    private static Set<String> possibleFinalStrings = new HashSet<String>();

    public VcfProcessor(List<VariantContext> leftVcfList) {
        this.leftVcfList = new ArrayList<VariantContext>(leftVcfList);
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
                        fatherNucleotideStringBuffer.replace(variantContext.getStart() - fatherStringShift,
                                variantContext.getEnd() - fatherStringShift + 1,
                                variantContext.getAlleles().get(1).getBaseString());
                        fatherStringShift += variantContext.getAlleles().get(0).getBaseString().length()
                                - variantContext.getAlleles().get(1).getBaseString().length();
                    } else {
                        motherNucleotideStringBuffer.replace(variantContext.getStart() - motherStringShift,
                                variantContext.getEnd() - motherStringShift + 1,
                                variantContext.getAlleles().get(0).getBaseString());
                        motherStringShift += variantContext.getAlleles().get(1).getBaseString().length()
                                - variantContext.getAlleles().get(0).getBaseString().length();
                    }
                } else {
                    StringBuilder fatherNucleotideStringBufferForCopy = new StringBuilder(fatherNucleotideStringBuffer);
                    fatherNucleotideStringBufferForCopy.replace(variantContext.getStart() - fatherStringShift,
                            variantContext.getEnd() - fatherStringShift + 1,
                            variantContext.getAlleles().get(1).getBaseString());
                    int fatherStringShiftForCopy = fatherStringShift +
                            variantContext.getAlleles().get(0).getBaseString().length()
                            - variantContext.getAlleles().get(1).getBaseString().length();

                    VcfProcessor vcfProcessorWithoutCurrentMutation = new VcfProcessor(this.leftVcfList);
                    vcfProcessorWithoutCurrentMutation.process(motherNucleotideStringBuffer.toString(),
                            fatherNucleotideStringBufferForCopy.toString(), motherStringShift,
                            fatherStringShiftForCopy);

                    motherNucleotideStringBuffer.replace(variantContext.getStart() - motherStringShift,
                            variantContext.getEnd() - motherStringShift + 1,
                            variantContext.getAlleles().get(1).getBaseString());
                    motherStringShift += variantContext.getAlleles().get(0).getBaseString().length()
                            - variantContext.getAlleles().get(1).getBaseString().length();
                }
            } else {
                motherNucleotideStringBuffer.replace(variantContext.getStart() - motherStringShift,
                        variantContext.getEnd() - motherStringShift + 1,
                        variantContext.getAlleles().get(1).getBaseString());
                fatherNucleotideStringBuffer.replace(variantContext.getStart() - fatherStringShift,
                        variantContext.getEnd() - fatherStringShift + 1,
                        variantContext.getAlleles().get(1).getBaseString());
                motherStringShift += variantContext.getAlleles().get(0).getBaseString().length()
                        - variantContext.getAlleles().get(1).getBaseString().length();
                fatherStringShift += variantContext.getAlleles().get(0).getBaseString().length()
                        - variantContext.getAlleles().get(1).getBaseString().length();
            }
        }

        possibleFinalStrings.add(motherNucleotideStringBuffer.toString());
        possibleFinalStrings.add(fatherNucleotideStringBuffer.toString());
    }

    public Set<String> getPossibleFinalStrings() {
        return possibleFinalStrings;
    }
}
