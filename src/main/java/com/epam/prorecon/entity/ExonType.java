package com.epam.prorecon.entity;

/**
 * Created by user on 7/28/2016.
 */
public enum ExonType {
    FIVE_UTR,
    EXON,
    START_CODON,
    CDS,
    STOP_CODON,
    THREE_UTR;

    public static ExonType asExonType(String str) {
        for (ExonType exonType : ExonType.values()) {
            if (exonType.name().equalsIgnoreCase(str))
                return exonType;
        }
        return null;
    }
}
