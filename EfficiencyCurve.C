#include <vector>

//ROOT header
#include "TVectorD.h"
#include "TGraph.h"
#include "TAxis.h"

//custom headers
#include "GGM_Analysis.h"



Int_t EfficiencyCurve(TString dst_files="");
void plot_efficiency(Int_t index, vector<double> array_hv_ch, vector<double> array_eff_ch);
Int_t plot_data_from_dst_list(TList* dst_folder);


void setAxisTitle(TGraph* hist, TString Xlabel, TString Ylabel);

/*
* in case of compilation with g++ we have a main function
*/
# ifndef __CINT__
int main(int argc, char* argv[]){
	
	//the first argument is the configuration file path
  return EfficiencyCurve( TString(argv[1]) );

}
# endif

Int_t EfficiencyCurve(TString dst_files){
	gGGM_Debug = 2;
	TList *dst_file_list;
	TSystemFile file_name = TSystemFile(dst_files.Data(), gSystem->DirName(dst_files.Data()));
	
	
	GGM_analysis_log( Form("\n\n/**********\nPlot Efficiency Curve started\n") );
	
	if( TString( file_name.GetName() ) == "" ){
		dst_files = "./";
		file_name = TSystemFile(gSystem->WorkingDirectory(), gSystem->DirName(gSystem->WorkingDirectory()));
	}
	file_name.Print();
	
	TString full_path = gSystem->ConcatFileName( gSystem->DirName(dst_files), gSystem->BaseName(dst_files) );
	
	//if is a directory, get the files in a lista
	if( file_name.IsDirectory() ){
		GGM_analysis_log( Form("Taking dst files from directory %s\n", full_path.Data()) );
		dst_file_list = list_directory_files(dst_files, ".dst");
	}else {
		GGM_analysis_log( Form("Taking dst files listed in %s\n", full_path.Data()) );
		dst_file_list = list_from_textfile(dst_files);
	}
	
	if( plot_data_from_dst_list(dst_file_list) == -1){
		GGM_analysis_log( Form("WARN: check dst files list at %s\n\n", full_path.Data()) );
		return -1;
	}
	else{
		return 0;
	}

	
}

/***
* plot efficiency curve, HV vs efficiency, for every channel from a list of dst files
*
* @param TString, folder path where every dst will be used
**/
Int_t plot_data_from_dst_list(TList* dst_file_list){
	
	
	if( dst_file_list == NULL || dst_file_list->GetEntries() == 0 ){
		GGM_analysis_log( Form("WARN: no dst file in this list %s\n") );
		return -1;
	}
	
		TListIter *dst_files_iter = (TListIter*)dst_file_list->MakeIterator(); //iterator object for the list
         
		TString filename;
        TSystemFile *file;
		
		double timestamp,ch_number,c,d,e,f,efficiency,hv_eff,i,l;
		
		vector< vector<double> > array_hv_ch(AVAILABLE_CHANNELS); //high voltage vector of AVAILABLE_CHANNELS of vector
		vector< vector<double> > array_eff_ch(AVAILABLE_CHANNELS); //efficiency vector of AVAILABLE_CHANNELS of vector
		
		
		int ch=0; //channel number
         while ((file=(TSystemFile*)dst_files_iter->Next())) {
            
			filename = file->GetName();
			
			GGM_analysis_log( Form("extracting data from file %s ...\n ", gSystem->ConcatFileName(file->GetTitle(), filename)) );

			ifstream dst_file( gSystem->ConcatFileName(file->GetTitle(), filename), ifstream::in); //open file in only-read mode
			
			if( dst_file.good() ){
				string line;
				ch=0; //inizilize channel number each time a dst file is opened
				cout << "file: " << filename << endl;
				dst_file >> timestamp >> ch_number >> c >> d >> e >> f >> efficiency >> hv_eff >> i >> l;
				while ( !dst_file.eof() && ch < AVAILABLE_CHANNELS) {
					  dst_file >> timestamp >> ch_number >> c >> d >> e >> f >> efficiency >> hv_eff >> i >> l;
					  cout << Form("lettura riga %g  -   %g\n", ch_number, efficiency);
					  array_hv_ch[ch_number].push_back( hv_eff/100 ); //HV
					  array_eff_ch[ch_number].push_back( efficiency/100 ); //efficiency
					  
					  ch++;
					  
				}
			}else{
				GGM_analysis_log( Form("WARN: skipped file %s ...\n ", gSystem->ConcatFileName(file->GetTitle(), filename)) );
			}
			

			cout << "fine file " << filename << endl <<endl;
			dst_file.close();
			

            GGM_analysis_log( Form("extraction finished.\n\n") );
         }
		 
		
		
		TGraph *gr;
		TCanvas *c1;
		TVectorD x_axis;
		TVectorD y_axis;
		
		//loop through vector to draw graph
		for( ch=0; ch < AVAILABLE_CHANNELS; ch++){
			
				
				if( array_hv_ch[ch].size() != 0){ //draw graph just for filled vector, exclude channel not in dst files
					
					plot_efficiency(ch+1, array_hv_ch[ch], array_eff_ch[ch]);

				}

				
		}
		
		//print all the canvases in PDF
		generate_images((TList*)gROOT->GetListOfCanvases(), "efficiency_curve.pdf");

		return 1;
}




void plot_efficiency(Int_t index, vector<double> vector_x, vector<double> vector_y){
	
	TGraph *gr;
	TCanvas *c1;
	TVectorD x_axis;
	TVectorD y_axis;
	
	cout << "size: " <<vector_x.size() <<endl;
	
	if( vector_x.size() == 0 || vector_y.size() == 0 ){
		GGM_analysis_log("WARN: no data for thsi channel.\n\n");
		return;
	}
	
	c1 = new TCanvas( Form("ch%d_graph", index), Form("ch%d graph", index), 800,800);
	
	x_axis.Use( vector_x.size(), &(vector_x[0]) ); //convert from C++ vector to TVectorD for TGraph compatibility
	y_axis.Use( vector_y.size(), &(vector_y[0]) );

	gr = new TGraph(x_axis, y_axis);
	gr->SetTitle(Form("ch%d graph", index));
	setAxisTitle(gr, "HV eff", "efficiency");
	
	gr->Draw("ACL*"); //draw graph in current pad
	
	//delete c1;
	
	return;
}


/**
*
* take an graph and set label on his axis
*
* @param TH1*
* @param TString label for abscissa
* @param TString label for ordinate
*
**/
void setAxisTitle(TGraph* graph, TString Xlabel, TString Ylabel){
	
	TAxis* asseX = graph->GetXaxis();
	asseX->SetTitle(Xlabel);
	asseX->SetTitleOffset(1.4);
	asseX->SetLabelSize(0.030);
	
	TAxis* asseY = graph->GetYaxis();
	asseY->SetTitle(Ylabel);
	asseY->SetTitleOffset(1.4);
	asseY->SetLabelSize(0.030);
	
	return;
}