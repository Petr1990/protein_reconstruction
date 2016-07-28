package com.epam.prorecon.entity;

public class ExonPosition {
    private String chromosome;
    private ExonType exonType;
    private int startIndex;
    private int endIndex;

    public ExonPosition(String chromosome, ExonType exonType, int startIndex, int endIndex) {
        this.chromosome = chromosome;
        this.exonType = exonType;
        this.startIndex = startIndex;
        this.endIndex = endIndex;
    }
}
