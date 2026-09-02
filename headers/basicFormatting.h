#ifndef BASICFORMATTING_H
#define BASICFORMATTING_H

#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TLatex.h>
#include <TPad.h>

//---Histogram formatting functions

void basicCanvasFormatting(TCanvas* c){
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);
}

void basicPaddedCanvasFormatting(TCanvas* c, TPad* pad1, TPad* pad2){
    // Canvas
    c->SetFillColor(0);
    c->SetFrameFillColor(0);
    c->SetTickx(1);
    c->SetTicky(1);

    // Top pad
    pad1->SetLeftMargin(0.098);
    pad1->SetRightMargin(0.036);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.02);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->SetFillColor(0);
    pad1->SetFrameFillColor(0);
    pad1->SetFrameLineWidth(1);

    // Bottom pad
    pad2->SetLeftMargin(0.098);
    pad2->SetRightMargin(0.036);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.30);

    pad2->SetTickx(1);
    pad2->SetTicky(0);
    pad2->SetFillColor(0);
    pad2->SetFrameFillColor(0);
    pad2->SetFrameLineWidth(1);
}

void basicHistFormatting(TH1D* hist){
    hist->GetXaxis()->CenterTitle(false);
    hist->GetYaxis()->CenterTitle(false);
    hist->GetXaxis()->SetTitleOffset(.9);
    hist->GetYaxis()->SetTitleOffset(1.);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetTitleSize(0.055);
}

void basicPaddedHistFormatting(TH1D* hist, bool isRatio = false){

    hist->SetTitle("");
    hist->SetStats(0);

    hist->GetXaxis()->CenterTitle(false);
    hist->GetYaxis()->CenterTitle(false);

    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);

    if (!isRatio) {

        // Main plot
        hist->GetXaxis()->SetTitleSize(0.055);
        hist->GetYaxis()->SetTitleSize(0.046);

        hist->GetXaxis()->SetLabelSize(0.0);  // Hide x labels
        hist->GetXaxis()->SetTitleSize(0.0);  // Hide x title

        hist->GetYaxis()->SetLabelSize(0.036);

        hist->GetYaxis()->SetTitleOffset(0.82);

    } else {

        // Ratio plot
        hist->GetYaxis()->CenterTitle(true);
        hist->GetXaxis()->SetTitleSize(0.1);
        hist->GetYaxis()->SetTitleSize(0.098);

        hist->GetXaxis()->SetLabelSize(0.08);
        hist->GetYaxis()->SetLabelSize(0.08);

        hist->GetXaxis()->SetTitleOffset(1.1);
        hist->GetYaxis()->SetTitleOffset(0.4);

        hist->GetXaxis()->SetTickLength(0.03);
        hist->GetYaxis()->SetTickLength(0.025);
    }

    hist->SetMarkerStyle(20);
}

void basicLegendFormatting(TLegend* leg){
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.044);
    leg->SetMargin(0.2);
    leg->SetEntrySeparation(0.04);
}

void drawLatexText(TString latexText, double x, double y, double TextSize){
// We can add any pT or eta selection. If no text is passed to the function the CMS Header will be drawn.
TLatex latex;
latex.SetNDC();              // For normalized coordinates
latex.SetTextSize(TextSize);
latex.SetTextFont(42);       // Helvetica
latex.SetTextAlign(11);      // Left-top aligned.
latex.DrawLatex(x, y, latexText);
}

#endif //BASICFORMATTING_H