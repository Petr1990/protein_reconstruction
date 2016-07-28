package com.epam.prorecon.entity;

public class WorkflowResult {
    private String appliedMutationsDnaString;
    private String mrnaString;
    private String mrnaWithoutIntronsString;
    private String proteinString;

    public String getAppliedMutationsDnaString() {
        return appliedMutationsDnaString;
    }

    public void setAppliedMutationsDnaString(String appliedMutationsDnaString) {
        this.appliedMutationsDnaString = appliedMutationsDnaString;
    }

    public String getMrnaString() {
        return mrnaString;
    }

    public void setMrnaString(String mrnaString) {
        this.mrnaString = mrnaString;
    }

    public String getMrnaWithoutIntronsString() {
        return mrnaWithoutIntronsString;
    }

    public void setMrnaWithoutIntronsString(String mrnaWithoutIntronsString) {
        this.mrnaWithoutIntronsString = mrnaWithoutIntronsString;
    }

    public String getProteinString() {
        return proteinString;
    }

    public void setProteinString(String proteinString) {
        this.proteinString = proteinString;
    }
}
