//c++ -O3 nid_pmt_waveforms_updated.cpp -lstdc++fs -o nid_analysis_updated.out `root-config --glibs --cflags`
// ./nid_analysis_updated.out <folder_with_files_to_analyze> <output_root_file> <number_of_rebinning_bins (typical 200)> <theshold_n_bins_above (typical -10)>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <experimental/filesystem>
#include <string>
#include <random>
#include <unistd.h>
#include <map>
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TLine.h"
#include "TROOT.h"
#include "TTree.h"
#include "lecroyparser.hpp"


namespace fs = std::experimental::filesystem;
using namespace std;

void print_histogram(TH1F *histo, string title, string x_axis, string y_axis, bool second_color = false, int sec_bin_min = 0, int sec_bin_max = 0);
void print_graph_lines (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax, TLine *l1, TLine *l2);
void print_graph_simple (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax);

tuple < double, double > calculate_RMS (int numb_points, shared_ptr<vector<double>> vec_ampl, int beginning);
void create_and_print_wf_graph_lines (string filename, shared_ptr<vector<double>> time, shared_ptr<vector<double>> ampl, double start, double end) ;
void create_and_print_wf_graph_simple (string filename, shared_ptr<vector<double>> time, shared_ptr<vector<double>> ampl, string tag);
vector <string> trim_name( string full_name, char delimiter);
void read_file_time_vs_amplitude (string file, string& names, shared_ptr<vector<double>>& time, shared_ptr<vector<double>>& amplitude);

void inject_fake_peaks (shared_ptr<vector<double>> time, shared_ptr<vector<double>> ampl, double threshold, double injection_ratio);

bool remove_noisy_files (string name, shared_ptr<vector<double>> t_vec, shared_ptr<vector<double>> a_vec) {

    double int_thresh = 3.e-7; // density of photons should be higher than this.
    int int_limit = 150; // signal values below 150 counts.
    bool skip_file = false;
    
    cout << 
    "\n\n\t*** Checking for noisy files..." << 
    "\n\nInteference density [t_diff]: " << int_thresh << 
    "\nInterference limit [#t_difs]: " << int_limit << endl;
    cout << "\nChecking file: " << name << endl;
    
    int intereference_counts = 0;
    double diff = 0;
    
    for (int k = 1; k < t_vec->size(); k++) {

        diff = ( (*t_vec)[k] - (*t_vec)[k-1] );
        
        if (diff <= int_thresh) intereference_counts++;

    }  cout << "intereference counts: " << intereference_counts << endl;

    if ( intereference_counts > int_limit){

        cout << "**Interference file**. Named: " << name << endl;
    }   

    return skip_file;
}

void analyse_histogram (TH1F *histo, double h_t_high_low, map<string, double>& histMap) {

    // Parameters

    double hist_thresh_high = h_t_high_low;
    double d_l_high = 2; //need 2 consecutive signals below X mV to start signal
    double hist_thresh_low = h_t_high_low;  
    double d_l_low = 2; //need 2 consecutive signals at X mV to stop signal

    int density_counts_start, density_counts_end;
    double area;
    bool beginning, ending;

    double signal_begin, signal_end;
    double time_window;

    double bin_value;
    double bin_width;
    int bin_begin, bin_end;

    cout << "\n** Histogram analysis **" << endl;

    density_counts_start = 0, density_counts_end = 0;
    beginning = false, ending = false;
    signal_begin = 0, signal_end = 0;
    bin_begin = 0, bin_end = 0;
    area = 0;

    bin_width = histo->GetXaxis()->GetBinWidth(1);

    for (int bin = 1; bin <= histo->GetNbinsX(); bin++){

        //Cross-check bins info
        // cout << "bin: " << bin << ", with value: " << histo->GetBinContent(bin) << ", centered on: " << histo->GetXaxis()->GetBinCenter(bin) << " ms." << endl;

        bin_value = histo->GetBinContent(bin);

        if (bin_value < hist_thresh_high) density_counts_start++;
        else density_counts_start = 0;

        if (density_counts_start  == d_l_high && beginning == false) {

            beginning = true;
            bin_begin = bin - d_l_high + 1;
            signal_begin = histo->GetXaxis()->GetBinCenter( bin_begin ) - bin_width/2.;
            cout << "Signal_begin: " << signal_begin << " ms" << endl;                
        }         

        if ( beginning == true){

            if (bin_value >= hist_thresh_low) density_counts_end++;
            else density_counts_end = 0;

            if ( density_counts_end  == d_l_low && ending == false ) {
            
                ending = true;
                bin_end = bin - d_l_low + 1;
                signal_end = histo->GetXaxis()->GetBinCenter( bin_end ) - bin_width/2.;
                cout << "Signal_end: " << signal_end << " ms" <<endl;
            }       
        }
    }

    if ( beginning == false ) {

        cout << "** beginning not found **" << endl;
        signal_begin = 0.;
        signal_end = 0.;
    }
    time_window = signal_end - signal_begin; 
    cout << "Signal width: " << time_window << endl;

    for (int bin_area = bin_begin; bin_area < bin_end; bin_area++){
        
        area += histo->GetBinContent( bin_area ) * bin_width;
    }
    cout << "Total area: " << abs(area) << endl;

    histMap["Time_window"] = time_window;
    histMap["Area"] = abs(area);
    histMap["bin_start"] = bin_begin;
    histMap["bin_end"] = bin_end;
    histMap["signal_begin"] = signal_begin;
    histMap["signal_end"] = signal_end;
}


int main(int argc, char**argv) {

    bool batch_mode = true;

    TApplication *myapp=new TApplication("myapp",0,0);
    
    if (batch_mode) gROOT->SetBatch(1);              // to avoid prints

    // =====================
    // Inputs
    // =====================
    
    string folder = argv[ 1 ];
    cout << "\nfolder: " << folder << "\n" << endl;

    string outputfile = argv[ 2 ]; 
    TFile* file_root = new TFile(outputfile.c_str(),"recreate");

    // Number of bins used in the re-binned and level of thresh to start/end signal
    int nBins = atof(argv[3]);                      // typical = 200
    double thresh_n_bins_above = atof(argv[4]);     // typical = -10


    // =====================
    // Analysis variables
    // =====================
    
    string name_of_file; 
    
    //  Threshold - Fixed to 6*RMS
    double threshold;       
    
    // Baseline and RMS
    double baseline, root_mean_square;

    // File identification
    // string identifier_name = "pmt_zoom--4th--400V_cm--00000"; // if one wants to read a specific file
    string identifier_name = "pmt_zoom";

    // Scan through files
    string filename_full;
    int files_read = 0;


    // =====================
    // Histograms (limits no longer critical since now there is a tree output)
    // =====================
    TH1F *hWaveform_averaged = new TH1F("","", nBins*2, -5.,8.);
    TH1F *hArea = new TH1F("","",120,0,300);
    TH1F *hTime_window = new TH1F ( "" , "", 100, 0., 10.);         
    TH1F *hPeaksHeight = new TH1F("","", 300, 0, 30);
    TH1F *hNumberOfPeaks = new TH1F("","", 100, 0, 1000);
    TH1F *hPeak_Freq= new TH1F ( "" , "", 1e3, 1e2, 1e7);
    TH1F *hRMS = new TH1F ( "" , "", 200, 0., 1.);
    TH1F *h6RMS = new TH1F ( "" , "", 200, 0., 6.);



    // =====================
    // Histogram information
    // =====================
    map<string, double > histoInfo;
    histoInfo["Area"] = 0;
    histoInfo["Time_window"] = 0;
    histoInfo["bin_start"] = 0;
    histoInfo["bin_end"] = 0;
    histoInfo["signal_begin"] = 0;
    histoInfo["signal_end"] = 0;

    // =====================
    // Tree output
    // =====================

    TTree *tree = new TTree("NIDTree", "Histogram Information");
    // Define variables to hold the values
    double tree_area = 0;           tree->Branch("Area",        &tree_area,         "Area/D");
    double tree_time_window = 0;    tree->Branch("Time_window", &tree_time_window,  "Time_window/D");
    double tree_bin_start = 0;      tree->Branch("bin_start",   &tree_bin_start,    "bin_start/D");
    double tree_bin_end = 0;        tree->Branch("bin_end",     &tree_bin_end,      "bin_end/D");
    double tree_signal_begin = 0;   tree->Branch("signal_begin",&tree_signal_begin, "signal_begin/D");
    double tree_signal_end = 0;     tree->Branch("signal_end",  &tree_signal_end,   "signal_end/D");
    double tree_baseline = 0;       tree->Branch("baseline",    &tree_baseline,     "baseline/D");
    double tree_6RMS = 0;           tree->Branch("RMS6",        &tree_6RMS,         "RMS6/D");
    double nPeaks = 0;              tree->Branch("nPeaks",      &nPeaks,            "nPeaks/D");

    // =====================
    // OPTIONS (main analysis controls)
    // =====================

    bool check_consecutive_points   = false;
    bool draw_raw_plot              = true;
    bool draw_lines_plot            = true;
    bool draw_above_noise_plot      = true;

    bool inject_peaks_for_testing   = false;       // Enable to test algorithm robustness
    double injection_ratio          = 0.1;         // 10% additional fake peaks in signal window

    // Number of files to read and plot
    int number_of_files_to_read = 6;    // if one doesn't want to read all the files. Offset to huge number (1e5) to read all files
    int number_of_files_to_plot = 6;    // number of files to make plots for
    int plots = 0;  

    int nRmsBins = 2000;    


    // =====================
    // Analysis start
    // =====================

    for (const auto & entry : fs::directory_iterator(folder)){
        filename_full = entry.path();
        if (filename_full.find(identifier_name) != string::npos) {

            if (entry.path().extension() != ".trc") continue;
            
            if (files_read >= number_of_files_to_read) break;

            // Retrieving data from file //
            shared_ptr<vector<double>> vec_time_raw = make_shared<std::vector<double>>();
            shared_ptr<vector<double>> vec_ampl_raw = make_shared<std::vector<double>>();

            read_file_time_vs_amplitude(filename_full, name_of_file, vec_time_raw, vec_ampl_raw);
            cout << "\n\tAnalysing file: " << name_of_file << " ...\n" << endl;

            // Other options // 
            shared_ptr<vector<double>> vec_time_noise = make_shared<std::vector<double>>();
            shared_ptr<vector<double>> vec_ampl_noise = make_shared<std::vector<double>>();

            // Calculating RMS noise//
            tie ( baseline, root_mean_square) = calculate_RMS( nRmsBins, vec_ampl_raw, 1);
            double RMS_6 = ( (  root_mean_square * -1. ) * 6. + baseline );
            hRMS->Fill(root_mean_square *1.e3);
            h6RMS->Fill(root_mean_square *6. *1.e3);
            cout << "** Basic waveform properties **" << endl;
            cout << "Baseline: " << baseline << "\nRms:" << root_mean_square << endl;
            cout << "6*RMS: " << RMS_6 << endl;
            threshold = RMS_6;

            // Starting histogram analysis // 
            TH1F *hWaveform = new TH1F("","", nBins, -5.,8.);

            double diff, freq;
            int diff_point;
            int points_above = 0;

            for (int k = 0; k < vec_time_raw->size(); k++){

                // get only the peaks above the noise //
                if ( (*vec_ampl_raw)[k] < threshold) {

                    // Cut on consecutive points above threshold //
                    // Avoid double-counting consecutive points above threshold
                    if ( check_consecutive_points == true) {
                    
                        if ( k > 0 && (*vec_ampl_raw)[k-1] < threshold ) continue;
                    }

                    hWaveform -> Fill ( (*vec_time_raw)[k] * 1.e3, (*vec_ampl_raw)[k] * 1.e3);
                    hPeaksHeight -> Fill ((*vec_ampl_raw)[k]*1.e3*(-1.));
                    hWaveform_averaged -> Fill ( (*vec_time_raw)[k] * 1.e3, (*vec_ampl_raw)[k] * 1.e3);

                    // To calculate the peak frequency // 
                    if (points_above > 0) {

                        diff = ( ((*vec_time_raw)[k]) - ((*vec_time_raw)[diff_point]) );
                        freq = 1. / diff;
                        hPeak_Freq -> Fill (freq);
                    }
                    
                    diff_point = k;
                    points_above++;
                }
                
                // To make the noise cut plot //
                if ( draw_above_noise_plot == true && files_read < number_of_files_to_plot) {
                    
                    if ( (*vec_ampl_raw)[k] < threshold) {

                        vec_time_noise->push_back((*vec_time_raw)[k]);
                        vec_ampl_noise->push_back((*vec_ampl_raw)[k]);
                    } 
                    else {
                        
                        vec_time_noise->push_back((*vec_time_raw)[k]);
                        vec_ampl_noise->push_back(0.);
                    }
                }
            }

            // Running the histogram analysis // 
            analyse_histogram(hWaveform, thresh_n_bins_above, histoInfo);

            hArea->Fill(histoInfo["Area"]);
            hTime_window->Fill(histoInfo["Time_window"]);
            hNumberOfPeaks->Fill(points_above);

            // Filling the tree //
            tree_area           = histoInfo["Area"];
            tree_time_window    = histoInfo["Time_window"];
            tree_bin_start      = histoInfo["bin_start"];
            tree_bin_end        = histoInfo["bin_end"];
            tree_signal_begin   = histoInfo["signal_begin"];
            tree_signal_end     = histoInfo["signal_end"];
            tree_baseline       = baseline;
            tree_6RMS           = RMS_6;
            nPeaks              = points_above;
            tree->Fill();

// ---------------------------- fake peak injection test starts ---------------------------- //

            // Optional: test with injected fake peaks
            if (inject_peaks_for_testing == true) {
                
                // Save original results
                map<string, double> orig_histoInfo = histoInfo;
                int orig_points_above = points_above;
                
                cout << "\n*** Testing with injected fake peaks ***" << endl;
                cout << "Original peak count: " << points_above << endl;
                
                // Make copy of noise vectors for testing
                shared_ptr<vector<double>> vec_time_test = make_shared<std::vector<double>>(*vec_time_noise);
                shared_ptr<vector<double>> vec_ampl_test = make_shared<std::vector<double>>(*vec_ampl_noise);
                
                // Inject fake peaks by replacing zeros
                inject_fake_peaks(vec_time_test, vec_ampl_test, threshold, injection_ratio);
                
                // Clear and re-analyze
                histoInfo.clear();
                hWaveform = new TH1F("","_modified", nBins, -5.,8.);
                points_above = 0;
                diff_point = 0;
                
                // Re-process with modified data
                for (int k = 0; k < vec_time_test->size(); k++) {
                    if ((*vec_ampl_test)[k] != 0.0) {
                        hWaveform -> Fill ( (*vec_time_test)[k] * 1.e3, (*vec_ampl_test)[k] * 1.e3);
                        if (points_above > 0) {
                            diff = ( ((*vec_time_test)[k]) - ((*vec_time_test)[diff_point]) );
                            if (diff > 1e-10) freq = 1. / diff;
                            hPeak_Freq -> Fill (freq);
                        }
                        diff_point = k;
                        points_above++;
                    }
                }
                
                // Re-analyze histogram
                analyse_histogram(hWaveform, thresh_n_bins_above, histoInfo);
                
                // Compare results
                cout << "\nComparison (Original vs With Injected Peaks):" << endl;
                cout << "  Area: " << orig_histoInfo["Area"] << " -> " << histoInfo["Area"] << endl;
                cout << "  Time_window: " << orig_histoInfo["Time_window"] << " -> " << histoInfo["Time_window"] << endl;
                cout << "  Peaks: " << orig_points_above << " -> " << points_above << endl;
                
                // Plot comparisons if we're plotting files
                if (plots < number_of_files_to_plot) {
                    create_and_print_wf_graph_simple(name_of_file + "_original", vec_time_noise, vec_ampl_noise, "peaks_original");
                    create_and_print_wf_graph_simple(name_of_file + "_injected", vec_time_test, vec_ampl_test, "peaks_injected");
                    
                    TH1F *hWaveform_orig = new TH1F("","_original", nBins, -5.,8.);
                    for (int k = 0; k < vec_time_noise->size(); k++) {
                        if ((*vec_ampl_noise)[k] != 0.0) {
                            hWaveform_orig -> Fill ( (*vec_time_noise)[k] * 1.e3, (*vec_ampl_noise)[k] * 1.e3);
                        }
                    }
                    
                    print_histogram(hWaveform_orig, (name_of_file + "_binned_original").c_str(), "t [ms]", "Amplitude [mV]", true, orig_histoInfo["bin_start"], orig_histoInfo["bin_end"]);
                    print_histogram(hWaveform, (name_of_file + "_binned_injected").c_str(), "t [ms]", "Amplitude [mV]", true, histoInfo["bin_start"], histoInfo["bin_end"]);
                    
                    delete hWaveform_orig;
                    plots++;
                }
                
                // delete hWaveform;
                histoInfo = orig_histoInfo;  // restore original for file output
            }

// ---------------------------- fake peak injection test ends ---------------------------- //

            // Doing some plots //
            if ( plots < number_of_files_to_plot && !inject_peaks_for_testing ) {

                print_histogram( hWaveform, (name_of_file+"_binned").c_str(),"t [ms]", "Amplitude [mV]", true, histoInfo["bin_start"], histoInfo["bin_end"]);

                if (draw_raw_plot == true) create_and_print_wf_graph_simple(name_of_file,vec_time_raw,vec_ampl_raw, "raw");
                if (draw_above_noise_plot == true) create_and_print_wf_graph_simple(name_of_file,vec_time_noise,vec_ampl_noise, "noise");
                if (draw_lines_plot == true) create_and_print_wf_graph_lines(name_of_file,vec_time_raw,vec_ampl_raw, histoInfo["signal_begin"], histoInfo["signal_end"]);
                
                plots++;
            }

            delete hWaveform;
            histoInfo.clear();
            files_read++;
        }
    }

    hWaveform_averaged -> Scale(-1.0 / files_read);
    print_histogram( hWaveform_averaged, "Average_Waveform","t [ms]", "Amplitude [mV]");

    print_histogram( hNumberOfPeaks,    "Number_of_peaks_per_waveform", "Number of peaks [#]",    "Counts [#]");
    print_histogram( hPeaksHeight,      "All_peaks_heights",            "Amplitude [mV]",         "Counts [#]");
    print_histogram( hPeak_Freq,        "All_peaks_frequencies",        "1/(#Deltat_peaks) [Hz]", "Counts [#]");

    print_histogram( hArea,         "Area",         "t [ms]",               "Amplitude [mV]");
    print_histogram( hTime_window,  "Time_window",  "t [ms]",               "Amplitude [mV]");
    print_histogram( hRMS,          "RMS",          "RMS amplitude [mV]",   "Counts [#]");
    print_histogram( h6RMS,         "6*RMS",        "6*RMS amplitude [mV]", "Counts [#]");


    tree->Write();
    cout << "\n\n\t**Finished**" << endl;
    cout << "\n\tData written into: " << outputfile <<  endl;

    file_root->Close();
    if (batch_mode == 0) myapp->Run();
    return 0;
}


// --- Helper functions ---



// void print_histogram(TH1F *histo, string title, string x_axis, string y_axis, bool second_color = false, int sec_bin_min = 0, int sec_bin_max = 0, bool save = false) {
void print_histogram(TH1F *histo, string title, string x_axis, string y_axis, bool second_color, int sec_bin_min, int sec_bin_max) {

    TCanvas *c = new TCanvas("","", 800, 400);
    c->cd();
    histo->SetTitle(title.c_str());
    histo->SetName(title.c_str());
    histo->SetLineColor(kAzure-5);
    histo->SetLineWidth(1);
    histo->SetFillStyle(3003);
    histo->SetFillColor(kAzure+5);
    histo->GetYaxis()->SetTitle(y_axis.c_str());
    histo->GetYaxis()->SetTitleOffset(1.0);
    histo->GetYaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetTitleOffset(1.0);
    histo->GetXaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetTitle(x_axis.c_str());
    histo->DrawCopy("hist");

    if ( second_color == true ){

        histo->GetXaxis()->SetRange(sec_bin_min,sec_bin_max);
        histo->SetLineColor(kRed);
        histo->SetFillColor(kRed);
        histo->DrawCopy("histsame");
        c->Write(title.c_str(),TObject::kWriteDelete);
    }

    if ( second_color == false ){
       
        histo->Write(title.c_str(),TObject::kWriteDelete);
    }
}

void print_graph_lines (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax, TLine *l1, TLine *l2){

    TCanvas *c = new TCanvas("","", 800, 400);
    c->cd();
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(x_axis.c_str());
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetXaxis()->SetTitleOffset(1);
    graph->GetYaxis()->SetTitle(y_axis.c_str());
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleOffset(1);
    graph->SetLineColor(kAzure-5);
    graph->SetMarkerColor(kAzure-5);
    graph->GetYaxis()->SetRangeUser(yMin,yMax);
    graph->Draw("apl");

    l1->SetLineColor(kGray+1);
    l1->SetLineWidth(2);
    l1->SetLineStyle(7);
    l1->Draw("same");
    l2->SetLineColor(kGray+1);
    l2->SetLineWidth(2);
    l2->SetLineStyle(7);                
    l2->Draw("same");

    c->SetName(title.c_str());
    c->Write(title.c_str(),TObject::kWriteDelete);
} 

void print_graph_simple (TGraph *graph, string title, string x_axis, string y_axis, double yMin, double yMax){

    TCanvas *c = new TCanvas("","", 800, 400);
    c->cd();
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(x_axis.c_str());
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetXaxis()->SetTitleOffset(1);
    graph->GetYaxis()->SetTitle(y_axis.c_str());
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleOffset(1);
    graph->SetLineColor(kAzure-5);
    graph->SetMarkerColor(kAzure-5);
    graph->GetYaxis()->SetRangeUser(yMin,yMax);
    graph->Draw("apl");
} 

tuple < double, double > calculate_RMS ( int numb_points, shared_ptr<vector<double>> vec_ampl, int beginning){

    double rms = 0;
    double bl = 0;

    for ( int i = 0 + beginning; i < numb_points + beginning; i ++) {

        bl += (*vec_ampl)[i] / numb_points;
    }

    for ( int i = 0 + beginning; i < numb_points + beginning; i ++) {

        rms += pow((*vec_ampl)[i] - bl,2) / numb_points;
    }

    rms = sqrt( rms );
    return make_tuple ( bl, rms );
}

void inject_fake_peaks (shared_ptr<vector<double>> time, shared_ptr<vector<double>> ampl, double threshold, double injection_ratio) {
    
    if (ampl->empty()) return;
    
    // Count existing peaks (non-zero amplitudes)
    int peak_count = 0;
    for (int i = 0; i < ampl->size(); i++) {
        if ((*ampl)[i] != 0.0) peak_count++;
    }
    
    if (peak_count == 0) return;
    
    // Calculate number of fake peaks to inject (10% of existing peaks)
    int num_to_inject = (int)(peak_count * injection_ratio);
    if (num_to_inject == 0) num_to_inject = 1;
    
    // Random number generator
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> dis(0, ampl->size() - 1);
    std::uniform_real_distribution<> amp_var(0.7, 1.3);
    
    // Find an existing peak to use as template
    double template_peak = threshold * 0.8;  // default fake peak amplitude
    for (int i = 0; i < ampl->size(); i++) {
        if ((*ampl)[i] != 0.0) {
            template_peak = (*ampl)[i];
            break;
        }
    }
    
    // Inject fake peaks by replacing zeros with peak amplitudes
    int injected = 0;
    int attempts = 0;
    while (injected < num_to_inject && attempts < num_to_inject * 10) {
        int random_idx = dis(gen);
        // Only replace zeros with fake peaks
        if ((*ampl)[random_idx] == 0.0) {
            (*ampl)[random_idx] = template_peak * amp_var(gen);
            injected++;
        }
        attempts++;
    }
    
    cout << "\n*** Injected " << injected << " fake peaks (" << (injection_ratio*100) << "%) replacing zeros ***" << endl;
}

void create_and_print_wf_graph_lines (string filename, shared_ptr<vector<double>> time, shared_ptr<vector<double>> ampl, double start, double end) {

    TGraph *gWaveform = new TGraph();
    string newname = filename + "_lines";

    TLine * line1 = new TLine(0.,-20.,0.,3.);
    line1->SetX1(start), line1->SetX2(start);
    TLine * line2 = new TLine(0.,-20.,0.,3.);
    line2->SetX1(end), line2->SetX2(end);

    for (int k = 0; k < time->size(); k++){

        gWaveform -> SetPoint ( k, (*time)[k] * 1.e3, (*ampl)[k] * 1.e3);
    }
    print_graph_lines(gWaveform, newname, "t [ms]", "Amplitude [mV]", -20.,3., line1, line2);
    
}

void create_and_print_wf_graph_simple (string filename, shared_ptr<vector<double>> time, shared_ptr<vector<double>> ampl, string tag) {

    TGraph *gWaveform = new TGraph();
    string newname = filename + "_" + tag;

    for (int k = 0; k < time->size(); k++){

        gWaveform -> SetPoint ( k, (*time)[k] * 1.e3, (*ampl)[k] * 1.e3);
    }
    print_graph_simple(gWaveform, newname, "t [ms]", "Amplitude [mV]", -20.,3.);
    gWaveform->SetName(newname.c_str());
    gWaveform->Write(newname.c_str(),TObject::kWriteDelete);
}

vector <string> trim_name( string full_name, char delimiter) {

    vector <string> tokens;
    stringstream check1(full_name);
    string intermediate;
    while(getline(check1, intermediate, delimiter)) tokens.push_back(intermediate);
    return tokens;
}

void read_file_time_vs_amplitude (string file, string& names, shared_ptr<vector<double>>& time, shared_ptr<vector<double>>& amplitude){

    string line, extension;
    ifstream myfile;
    double t, ampl;

    myfile.open(file.c_str());
    vector <string> name_trim = trim_name(file, '/');
    names = name_trim.back().c_str(); 

    if ( trim_name(name_trim.back(),'.')[1].compare("txt") == 0) {
        extension = "txt";
    }
    else if ( trim_name(name_trim.back(),'.')[1].compare("trc") == 0) {
        extension = "trc";
    }

    if ( extension.compare("txt") == 0 ) {

        for ( int j = 0; j < 5; j++)  getline( myfile, line );  //ignore first 5 lines 
        while ( getline ( myfile, line ) ) {

            istringstream iss ( line );	//creates string consisting of a line
            string token;

            getline (iss, token, '\t');
            t = (double) stod(token);
            time->push_back(t);

            getline (iss, token, '\t');
            ampl = (double) stod(token);
            amplitude->push_back(ampl);
        }
    }  

    else if ( extension.compare("trc") == 0 ) {   
        
        lecroyparser::lecroy_wavedesc_2_3 header;
        int16_t* waveform = nullptr;
        bool res = lecroyparser::read(file, header, waveform);

        double x_offset = header.horiz_offset;
        double x_interval = header.horiz_interval;
        double y_gain = header.vertical_gain;
        double y_offset = header.vertical_offset;

        for ( long int i = 0; i < header.wave_array_count; ++i) {
            
            time->push_back(x_offset + x_interval * i );
            amplitude->push_back( waveform[i] * y_gain - y_offset);
        }

        delete[] waveform;
    }

    myfile.close();
}