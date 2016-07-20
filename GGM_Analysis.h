/*****
*
* header file with common utility function
* between GGM Analysys software
*
*******/
#ifndef __GGM_Analysis__ //header guard lock
#define __GGM_Analysis__

#include <iostream>
#include <fstream>

using namespace std;


//ROOT headers
#include "TROOT.h"
#include "TEnv.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TError.h"





#define AVAILABLE_CHANNELS 16 //how many ADC channels
#define MINIMUM_ENTRIES 5000 //less then this number of events will not be analyzed


/*
* struct type definition for managing the configuration file parameters
*/
typedef struct{
   TString pedestal_filename;
   TString total_signal_filename;
   TString dst_filename;
   TString output_filename;
   TString excluded_channels;
   Int_t debug_mode;
}AnalysisConfig;



void parse_config_file(TString, AnalysisConfig*);
Int_t find_number(TString needle, TString haystack);
Int_t find_number(Int_t needle, TString haystack);

void GGM_analysis_log(TString message);
void GGM_analysis_log(const char* message);

void printArray(Double_t arr[], int size);

void generate_images(TList *, TString);

TList* list_directory_files(TString dir_path, TString ext="");
TList* list_from_textfile(TString textfile);

/******
*
* global variable
*******/
Int_t gGGM_Debug = 0; //global variable for debug option, 0 print nothing in output




/**
* loop through a list of canvases and print each one in a single PDF multipage
*
* @param TList*, puntatore a una lista, gli elementi devono essere istanze della classe TCanvas
**/
void generate_images(TList *canvas_list, TString file_path){
	
	gErrorIgnoreLevel=1; //impedisce che vengano stampati messaggi di avvisi sulla linea di comando
	
	if( canvas_list == NULL || canvas_list->GetEntries() == 0){
		GGM_analysis_log("No canvas to save\n\n");
		return;
	}
	
	TCanvas *canvas;
	TPaveText *t; 
	TString name = gSystem->BaseName(file_path.Data());;

	Int_t index=0;
	Int_t list_entries = canvas_list->GetEntries();

	
	//create directory path, if not exists
   TString dir = gSystem->DirName(file_path);
   if( gSystem->AccessPathName(dir) ) //retun true if NOT exists
      if( gSystem->mkdir( dir , kTRUE) == -1)	
         GGM_analysis_log( Form("ERROR: can't create directory %s: %s\n", dir.Data()) );

	
	canvas_list->At(0)->Print( Form("%s(", file_path.Data()) ); //first canvas
	
	TListIter *iter = (TListIter*)canvas_list->MakeIterator(); //loop over middle canvases
	for(index=1; index < list_entries - 1; index++) {
		canvas = (TCanvas*)canvas_list->At(index);
		canvas->cd(1);
		
		
		//add file name to canvas
		//t = new TPaveText(0.1,0.004,0.39,0.05,"nbNDC");
		//t->AddText( file_path.Data() );
		//t->Draw();
		  
		canvas->Print( Form("%s", file_path.Data()) );
	}
	
	canvas_list->At(list_entries-1)->Print( Form("%s)", file_path.Data()) ); //last canvas
	return;
}


/**
* @use global gGGM_Debug
* @param TString with message to display only if in debug mode
*
**/
void GGM_analysis_log(TString message){
	
	if( gGGM_Debug > 0 )
		cout << message;

}

/**
* second signature for GGM_analysis_log
* @use global gGGM_Debug
* @param char* array with message to display only if in debug mode
*
**/
void GGM_analysis_log(const char* message){
	
	GGM_analysis_log( TString(message) );

}


/**
*
* @use global gGGM_Debug
* @param TString* config file path
* @param AnalysisConfig*, pointer to the struct will be filled with configuration parameters
*
**/
void parse_config_file(TString file_path, AnalysisConfig* config){

	//reading config file
	TEnv analysis_config;
	if( !gSystem->AccessPathName(file_path.Data()) ) //if file exists, load it, else it will use default values
		analysis_config.ReadFile(file_path.Data(), kEnvLocal );

	if( gDebug > 0 ){
	  cout << Form("configuration file path: %s\n\n", gSystem->ConcatFileName(gSystem->pwd(), file_path.Data()));
	}
	
	// struct to save configuration parameters,
	// the second parameters are default values
	config->pedestal_filename = analysis_config.GetValue("pedestal-file", "temp1.raw");
	config->total_signal_filename = analysis_config.GetValue("total-signal-file", "temp2.raw");
	config->dst_filename = analysis_config.GetValue("dst-file", "temp.dst");
	config->output_filename = analysis_config.GetValue("output-file", "temp.pdf");
	config->excluded_channels = analysis_config.GetValue("excluded_channels", "13");
	config->debug_mode = analysis_config.GetValue("debug_mode", 0);
	
	if( config->debug_mode > 0){
		gGGM_Debug = config->debug_mode;
	}
	
}

/*
* @param needle stringa da cercare
* @param haystack stinga in cui cercare needle
* @return bool vero in caso needle è trovato in haystack
**/
Int_t find_number(TString needle, TString haystack){

	TString tok;
	Ssiz_t from = 0;
	while (haystack.Tokenize(tok, from, ",")) {
	   if( needle.EqualTo(tok) ){
         return 1;
      }
	}
	
	return 0;
}

/*
* @param needle numero intero da cercare
* @param haystack stinga in cui cercare needle
* @return bool vero in caso needle è trovato in haystack
**/
Int_t find_number(Int_t needle, TString haystack){
	
	return find_number( Form("%d",needle), haystack);
}


/****
* utility for printing array
*****/
void printArray(double arr[], int size) {
    for ( int i = 0; i < size; i++ ) {
        cout << arr[i] << ' ';
    }
    cout << endl;
}

/****
* utility for printing vector
*****/
void printArray(vector<double> arr, int size) {
    for ( int i = 0; i < size; i++ ) {
        cout << arr[i] << ' ';
    }
    cout << endl;
}



/**
* return a TList containing all files in the directory  dir_path with extension ext
*
* @param directory path
* @param extension to filter out
* @return a TList* to the list of file or NULL
**/
TList* list_directory_files(TString dir_path, TString ext){

   TSystemDirectory dir( dir_path, dir_path);
   TList *files = dir.GetListOfFiles();

   if (files) { //directory colud be empty
		TSystemFile *file;
		TString fname;
		TListIter *iter = (TListIter*)files->MakeIterator(); //iterator object for the list
	
		//loop throw the dir files to remove file without the specified extension
		while ((file=(TSystemFile*)iter->Next())) {
		 fname = file->GetName();
		 if (file->IsDirectory() || !fname.EndsWith(ext)) {
			//cout << fname.Data() << endl;
			files->Remove(file);
		 }
	   }
   }
   
   return files;
}

/**
* return a TList containing all files in listed in the text file
* each line is supposed to be a file path to a dst file 
*
* @param text file path
* @return a TList* to the list of file or NULL
**/
TList* list_from_textfile(TString textfile){
	
	TList* list_of_files = new TList();
	TSystemFile *file;
	
	ifstream textfile_stream( textfile.Data(), ifstream::in); //open file in only-read mode	
	TString file_path;
	
	if( !textfile_stream ){
		GGM_analysis_log( Form("ERROR: Can't read file %s\n\n", textfile.Data()) );
		return NULL;
	}
	
	while ( !textfile_stream.eof() ) {
		
		textfile_stream >> file_path;
		file = new TSystemFile(file_path.Data(),"");
		if( !file->IsDirectory() && !gSystem->AccessPathName(file_path.Data()) )
			list_of_files->Add(file);

	}

	textfile_stream.close();
	
	return list_of_files;
}

#endif