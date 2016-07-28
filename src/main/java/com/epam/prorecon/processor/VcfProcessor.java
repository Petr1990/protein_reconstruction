package com.epam.prorecon.processor;

import com.epam.prorecon.entity.ExonPosition;
import com.epam.prorecon.entity.ExonType;
import com.epam.prorecon.entity.WorkflowResult;
import htsjdk.variant.variantcontext.VariantContext;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

public class VcfProcessor {
    private static final Logger LOGGER = LoggerFactory.getLogger(VcfProcessor.class);

    private static final Map<Character, Character> COMPLEMENT_MRNA_NUCLEOTIDE_MAP =
            Collections.unmodifiableMap(new HashMap<Character, Character>() {{
                put('A', 'U'); put('T', 'A'); put('C', 'G'); put('G', 'C');
            }});
    private static final Map<String, String> NUCLEOTIDE_TRIPLET_AMINO_ACID_MAP =
            Collections.unmodifiableMap(new HashMap<String, String>() {{
                put("UUU", "F"); put("UUC", "F"); put("UUA", "L"); put("UUG", "L"); put("CUU", "L"); put("CUC", "L");
                put("CUA", "L"); put("CUG", "L"); put("AUU", "I"); put("AUC", "I"); put("AUA", "I"); put("AUG", "M");
                put("GUU", "V"); put("GUC", "V"); put("GUA", "V"); put("GUG", "V"); put("UCU", "S"); put("UCC", "S");
                put("UCA", "S"); put("UCG", "S"); put("CCU", "P"); put("CCC", "P"); put("CCA", "P"); put("CCG", "P");
                put("ACU", "T"); put("ACC", "T"); put("ACA", "T"); put("ACG", "T"); put("GCU", "A"); put("GCC", "A");
                put("GCA", "A"); put("GCG", "A"); put("UAU", "Y"); put("UAC", "Y"); put("UAA", "$"); put("UAG", "$");
                put("CAU", "H"); put("CAC", "H"); put("CAA", "Q"); put("CAG", "Q"); put("AAU", "N"); put("AAC", "N");
                put("AAA", "K"); put("AAG", "K"); put("GAU", "D"); put("GAC", "D"); put("GAA", "E"); put("GAG", "E");
                put("UGU", "C"); put("UGC", "C"); put("UGA", "$"); put("UGG", "W"); put("CGU", "R"); put("CGC", "R");
                put("CGA", "R"); put("CGG", "R"); put("AGU", "S"); put("AGC", "S"); put("AGA", "R"); put("AGG", "R");
                put("GGU", "G"); put("GGC", "G"); put("GGA", "G"); put("GGG", "G");
            }});

    private List<VariantContext> leftVcfList;
    private static Set<String> possibleFinalStrings = new HashSet<>();
    private static Collection<WorkflowResult> workflowResults = new ArrayList<>();

    public VcfProcessor(List<VariantContext> leftVcfList) {
        this.leftVcfList = new ArrayList<>(leftVcfList);
    }

    public void process(String motherNucleotideString, String fatherNucleotideString, int motherStringShift,
                          int fatherStringShift, List<ExonPosition> motherExonPositions,
                        List<ExonPosition> fatherExonPositions)
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
                                fatherStringShift, fatherExonPositions);
                    } else {
                        motherStringShift = applyMutation(variantContext, motherNucleotideStringBuffer,
                                motherStringShift, motherExonPositions);
                    }
                } else {
                    StringBuilder fatherNucleotideStringBufferForCopy = new StringBuilder(fatherNucleotideStringBuffer);
                    int fatherStringShiftForCopy = applyMutation(variantContext, fatherNucleotideStringBufferForCopy,
                            fatherStringShift, fatherExonPositions);

                    VcfProcessor vcfProcessorWithoutCurrentMutation = new VcfProcessor(leftVcfList);
                    vcfProcessorWithoutCurrentMutation.process(motherNucleotideStringBuffer.toString(),
                            fatherNucleotideStringBufferForCopy.toString(), motherStringShift,
                            fatherStringShiftForCopy, new ArrayList<>(motherExonPositions),
                            new ArrayList<>(fatherExonPositions));

                    motherStringShift = applyMutation(variantContext, motherNucleotideStringBuffer, motherStringShift,
                            motherExonPositions);
                }
            } else {
                motherStringShift = applyMutation(variantContext, motherNucleotideStringBuffer, motherStringShift,
                        motherExonPositions);
                fatherStringShift = applyMutation(variantContext, fatherNucleotideStringBuffer, fatherStringShift,
                        fatherExonPositions);
            }
        }

        addNewResult(motherNucleotideStringBuffer, motherExonPositions, motherStringShift);
        addNewResult(fatherNucleotideStringBuffer, fatherExonPositions, fatherStringShift);
    }

    public void addReferenceResult(String nucleotideString, List<ExonPosition> exonPositions, int shift) {
        addNewResult(new StringBuilder(nucleotideString), exonPositions, shift);
    }

    private void addNewResult(StringBuilder nucleotideStringBuffer, List<ExonPosition> exonPositions, int shift) {
        if (!possibleFinalStrings.contains(nucleotideStringBuffer.toString())) {
            possibleFinalStrings.add(nucleotideStringBuffer.toString());

            Collections.sort(exonPositions, (ep1, ep2) -> ep1.getStartIndex() - ep2.getStartIndex());
            String nonInvertedMrnaString = generateNonInvertedMrnaString(nucleotideStringBuffer.toString());
            ExonPosition relativeStartCodonPosition = new ExonPosition();
            String invertedMrnaWithoutIntronsString = generateInvertedMrnaString(
                    new StringBuilder(generateMrnaWithoutIntronsString(nonInvertedMrnaString, exonPositions,
                            relativeStartCodonPosition, shift)));

            WorkflowResult workflowResult = new WorkflowResult();
            workflowResult.setAppliedMutationsDnaString(nucleotideStringBuffer.toString());
            workflowResult.setMrnaString(generateInvertedMrnaString(new StringBuilder(nonInvertedMrnaString)));
            workflowResult.setMrnaWithoutIntronsString(invertedMrnaWithoutIntronsString);
            workflowResult.setProteinString(generateProteinString(invertedMrnaWithoutIntronsString,
                    relativeStartCodonPosition.getStartIndex()));

            workflowResults.add(workflowResult);
        }
    }

    private int applyMutation(VariantContext variantContext, StringBuilder nucleotideStringBuilder, int shift,
                              List<ExonPosition> exonPositions) {
        nucleotideStringBuilder.replace(variantContext.getStart() - shift, variantContext.getEnd() - shift + 1,
                variantContext.getAlleles().get(1).getBaseString());

        int mutationShift = variantContext.getAlleles().get(0).getBaseString().length()
                - variantContext.getAlleles().get(1).getBaseString().length();
        int mutationShiftStart = variantContext.getStart() - shift +
                variantContext.getAlleles().get(1).getBaseString().length();
        int mutationShiftEnd = variantContext.getEnd() - shift + 1;

        for (ExonPosition exonPosition : exonPositions) {
            if (!(mutationShiftEnd < exonPosition.getStartIndex() || mutationShiftStart > exonPosition.getEndIndex())) {
                exonPosition.setStartIndex(Math.min(mutationShiftStart, exonPosition.getStartIndex()));
                exonPosition.setEndIndex(Math.max(exonPosition.getEndIndex() - mutationShiftEnd
                        + mutationShiftStart, mutationShiftStart));
            }
        }

        return shift + mutationShift;
    }

    private String generateNonInvertedMrnaString(String nucleotideString) {
        StringBuilder mrnaStringBuilder = new StringBuilder();
        
        for (int i = 0; i < nucleotideString.length(); i++){
            char nucleotide = nucleotideString.charAt(i);
            Character complementMrnaNucleotide = COMPLEMENT_MRNA_NUCLEOTIDE_MAP.get(nucleotide);
            if (complementMrnaNucleotide != null) {
                mrnaStringBuilder.append(complementMrnaNucleotide);
            } else {
                LOGGER.error("Failed to find complement mrna nucleotide to " + nucleotide);
            }
        }

        return mrnaStringBuilder.toString();
    }

    private String generateInvertedMrnaString(StringBuilder nonInvertedMrnaStringBuilder) {
        return nonInvertedMrnaStringBuilder.reverse().toString();
    }

    private String generateMrnaWithoutIntronsString(String nonInvertedMrnaString, List<ExonPosition> exonPositions,
                                                    ExonPosition relativeStartCodonPosition, int shift) {
        StringBuilder nucleotideStringBuilder = new StringBuilder();

        ExonPosition startCodon = new ExonPosition();
        ExonPosition exonWithStartCodon = new ExonPosition();
        Optional<ExonPosition> optionalStartCodon = exonPositions.stream().filter(exonPosition ->
                ExonType.START_CODON.equals(exonPosition.getExonType())).findFirst();
        if (optionalStartCodon.isPresent()) {
            startCodon.setExonPosition(optionalStartCodon.get());
            Optional<ExonPosition> optionalExonWithStartCodon = exonPositions.stream().filter(exonPosition ->
                    ExonType.EXON.equals(exonPosition.getExonType())
                            && exonPosition.getStartIndex() <= startCodon.getStartIndex()
                            && exonPosition.getEndIndex() >= startCodon.getEndIndex()).findFirst();
            if (optionalExonWithStartCodon.isPresent()) {
                exonWithStartCodon.setExonPosition(optionalExonWithStartCodon.get());
            }
        }

        exonPositions.stream().filter(exonPosition -> ExonType.EXON.equals(exonPosition.getExonType())).
                forEachOrdered(exonPosition -> {
                    if (exonWithStartCodon.equals(exonPosition)) {
                        relativeStartCodonPosition.setExonPosition(startCodon);
                        relativeStartCodonPosition.setStartIndex(nucleotideStringBuilder.length()
                                + startCodon.getStartIndex() - exonWithStartCodon.getStartIndex());
                        relativeStartCodonPosition.setEndIndex(nucleotideStringBuilder.length()
                                + startCodon.getEndIndex() - exonWithStartCodon.getStartIndex());
                    }

                    nucleotideStringBuilder.append(nonInvertedMrnaString.substring(exonPosition.getStartIndex() - shift,
                            exonPosition.getEndIndex() - shift));
                });

        return generateInvertedMrnaString(nucleotideStringBuilder);
    }

    private String generateProteinString(String invertedMrnaWithoutIntrons, int startCodonIndex) {
        StringBuilder proteinStringBuilder = new StringBuilder();

        if (!"AUG".equals(invertedMrnaWithoutIntrons.substring(startCodonIndex, startCodonIndex + 3))) {
            return "";
        }

        for (int i = startCodonIndex; i + 2 < invertedMrnaWithoutIntrons.length(); i += 3) {
            String aminoAcid = NUCLEOTIDE_TRIPLET_AMINO_ACID_MAP.get(invertedMrnaWithoutIntrons.substring(i, i + 3));
            if ("$".equals(aminoAcid)) {
                return proteinStringBuilder.toString();
            }
            proteinStringBuilder.append(aminoAcid);
        }

        return proteinStringBuilder.toString();
    }

    public Collection<WorkflowResult> getWorkflowResults() {
        return workflowResults;
    }
}
