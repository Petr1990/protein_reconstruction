package com.epam.prorecon.processor;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;

public abstract class AbstractProcessor {
    public abstract void process(String motherNucleotideString, String fatherNucleotideString, int motherStringShift,
                                   int fatherStringShift)
            throws CloneNotSupportedException;
}
