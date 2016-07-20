

//ROOT headers
#include "TH1F.h"
#include "TNtuple.h"
#include "TTimeStamp.h"
#include "TDatime.h"
#include "TBenchmark.h"
#include "TMath.h"

#include "TStopwatch.h"

//custom headers
#include "GGM_Analysis.h"

#define BIN_WIDTH 3.0 //fixed bin width

#define RMS_MIN 3 //minimum RMS to plot an histogram


void GGM_Analysis(TString config_filename = "ggm-analysis.conf");
Double_t analyze_channel_efficiency(Int_t ch_number, TNtuple*,TNtuple*);

Int_t check_valid_channel(Int_t, TNtuple*);
TH1F *RemoveOutliers(TNtuple* ntupla, Int_t ch_number);
Double_t FindFirstZeroBeforeMaximun(TH1* hist);
Double_t RoundUp(Double_t i, Double_t n);

TNtuple* populate_tree(TString file, TString title);
void print_extra_info(TH1 *hist);

Double_t efficiency_calc(TH1 *sgn_ped_diff, TH1 *sgn);
void update_dst_file(int index, double efficiency, TString filename = "");
void set_negative_to_zero(TH1* hist);
void set_negative_bin_to_zero(TH1* hist);

void setAxisTitle(TH1* hist, TString Xlabel, TString Ylabel);

double read_hv(TString filename, int k);


/*
* in case of compilation with g++ we have a main function
*/
# ifndef __CINT__
int main(int argc, char* argv[]){
	
	//the first argument is the configuration file path
  GGM_Analysis( TString(argv[1]) );
  return 0;
}
# endif


/**********
* @param TString path to config file
*
***********/
void GGM_Analysis(TString config_filename){
	
	TStopwatch execution_time;
	execution_time.Start();


	//read and parse configuration file
	AnalysisConfig config;
	parse_config_file(config_filename, &config);

	TNtuple *channels_sgn_tot;
	TNtuple *channels_sgn_ped;
	 
			
   //reading raw file for total signal
   channels_sgn_tot = populate_tree(config.total_signal_filename, "total_signal" );
   if( channels_sgn_tot == NULL ){
	   GGM_analysis_log( Form("Can't read file %s\nAnalysis aborted.\n\n", config.total_signal_filename.Data()) );
	   return;
   }

   //reading raw file for pedestal
   channels_sgn_ped = populate_tree(config.pedestal_filename, "pedestal" );
   if( channels_sgn_ped == NULL ){
	   GGM_analysis_log( Form("Can't read file %s\nAnalysis aborted.\n\n", config.pedestal_filename.Data()) );
	   return;
   }
   
   if( channels_sgn_ped->GetEntries() < MINIMUM_ENTRIES-(MINIMUM_ENTRIES*0.01) ){
	  GGM_analysis_log( Form("WARN: too few events in %s\nAnalysis aborted.\n", config.pedestal_filename.Data()) );
	  return;
   }
	   

	 
   Double_t efficiency = 0;
      

    //show the statistics box in histograms plot
   gStyle->SetOptStat("nemruo");


/****
* 2)
* loop over all AVAILABLE_CHANNELS of the ADC
*****/
for(Int_t i=1; i <= AVAILABLE_CHANNELS; i++){
	
	
	// if file is in this list, dont' analyze it
	if( !find_number(i, config.excluded_channels) ){
			
			// check for broken or power off channel
		   if( !check_valid_channel(i, channels_sgn_ped) ){
				GGM_analysis_log( Form("Channel %d not valid, skipped\n\n", i) );
				efficiency = 0; //salvo lo stesso l'efficenza del canale come valore zero
			  
		   }else{
				GGM_analysis_log( Form("Channel %d analysis...\n\n", i) );
				efficiency = analyze_channel_efficiency(i, channels_sgn_ped,channels_sgn_tot);
		   }

			//udpdate dst file with efficiency for this channel
			update_dst_file(i,efficiency, config.dst_filename);
	}else{
		GGM_analysis_log( Form("Channel %d excluded in configuration file, skipped\n\n", i) );
	}

}// end loop throw channels


	//loop throw all canvases and save a multipage PDF
	generate_images((TList*)gROOT->GetListOfCanvases(), config.output_filename);


	execution_time.Stop();
	execution_time.Print();
	

	return;
} //analysys end


/**
* check if the channel ch of the ntuple is valid:
*   1) check RMS without outliers to check id ADC channel was off or broken
*    2) check if data entries are enoguh to be analysed
*
* @param Int_t ch is the number of the channel
* @param Int_t ch is the Ntupla to get the channel
* @return boolean: if valid for analysis return true, else false 
**/
Int_t check_valid_channel(Int_t ch, TNtuple* ntuple){
   
   Double_t mean = 0;
   TH1F *ch_temp = NULL; //istogramma temporaneo per il piedistallo
	
	
	
   //return an histogram with a new range to remove outliers
   ch_temp = RemoveOutliers(ntuple, ch);
   if(ch_temp == NULL){
	   return kFALSE;
   }
   
   if( ch_temp->GetRMS() > RMS_MIN ){
	   return kTRUE;
   }else{
	   return kFALSE;
   }

}


/**
* esegue la procedura di sovrapposizione e fit dell'istogramma
* restituisce il valore dell'efficenza del canale
**/
Double_t analyze_channel_efficiency(Int_t ch_number, TNtuple* channels_sgn_ped, TNtuple* channels_sgn_tot){
	

	Int_t i = ch_number; //short name for channel index
   
   Double_t efficiency = 0; //variable to return
	
	
	/**
	* 1)
	* return a new histogram with optimal range and bin width 
	*/
   TH1F *ch_ped_temp = RemoveOutliers(channels_sgn_ped, i);
   


   /**
   * 3.1)
   * creo istogramma del piedistallo con i nuovi parametri partendo dalla TNtuple
   **/

      //piedistallo
   Double_t center = ch_ped_temp->GetBinCenter( ch_ped_temp->GetMaximumBin() ); //get the maximum, same as the mean for gaussian distibution
   
   //to fit the range and accomodate the total signal and the pedestal in the same range, calculate the lower limit as the most left value of the pedestal
   // if RMS is too big, cut it off
   Double_t new_rms = 15;
   if( ch_ped_temp->GetRMS() > 20.0 ){
	   new_rms = 7;
   }
   
	Double_t bin_width = BIN_WIDTH;
	Double_t lower_limit = center-5*ch_ped_temp->GetRMS();
	Double_t upper_limit = center+new_rms*ch_ped_temp->GetRMS();
		
	Int_t nbinsx = (upper_limit-lower_limit) / bin_width;
	
	lower_limit = RoundUp(lower_limit, bin_width);
	upper_limit = RoundUp(upper_limit, bin_width);

	
	/****
	* TEST
	* try to get optminal parameter directly from ch_ped_temp after RemoveOutliers
	*
	* Int_t nbinsx = ch_ped_temp->GetNbinsX()*2;
	* Double_t bin_width = ch_ped_temp->GetBinWidth(1);
	* Double_t lower_limit = ch_ped_temp->GetBinLowEdge(1);
	* Double_t upper_limit = (bin_width*nbinsx)+lower_limit;
	**/
  
	
	//creata canvas object to save histograms as images
    TCanvas *canvas = new TCanvas(Form("c%d_canvas", i),Form("canvas channel %d", i),800,800);
	canvas->Divide(1,1);
	canvas->cd(1);
   
	/**
	* 2)
	* pedestal histogram
	*  draw histogram of channel_%d minus zero from ntuple in a custom binning
	**/
   channels_sgn_ped->Draw(Form("(channel_%d-%g)>>ch%d_ped(%d,%g,%g)",i, center,i,nbinsx, lower_limit-center, upper_limit-center),""); 
   
	TH1F* ped_histogram = (TH1F*)gDirectory->FindObject( Form("ch%d_ped",i) );
	ped_histogram->SetLineColor(kRed);
	setAxisTitle(ped_histogram, "ADC charge", "count"); //set axis labels
	//print_extra_info(ped_histogram); //print stats infomation

   

	/**
	* 3+4)
	* draw total signal in an histogram with the binning as pedestal
	**/
	channels_sgn_tot->Draw(Form("(channel_%d-%g)>>ch%d_tot(%d,%g,%g)",i, center,i,nbinsx, lower_limit-center, upper_limit-center),"","same");

	TH1F* sgn_tot_histogram = (TH1F*)gDirectory->FindObject( Form("ch%d_tot",i) );
	sgn_tot_histogram->SetLineColor(kBlue);
   //print_extra_info(sgn_tot_histogram); //print stats infomation

   /**
   * 5)
   * scale the total signal
   **/
   Double_t scale_factor = ped_histogram->GetBinContent( ped_histogram->GetMaximumBin() ) / sgn_tot_histogram->GetBinContent( ped_histogram->GetMaximumBin() ) ;
   sgn_tot_histogram->Scale(scale_factor);



   /**
   * 6)
   * calculate the difference between total signal and pedestal
   **/
   TH1F *sgn_diff = (TH1F*)sgn_tot_histogram->Clone();
   sgn_diff->SetName(Form("ch%d_diff", i));
   sgn_diff->Add(ped_histogram,-1); //difference calculation

   set_negative_to_zero(sgn_diff); //set to zero negative abscissa because they have no meaning
   set_negative_bin_to_zero(sgn_diff);

   sgn_diff->SetTitle(Form("channel %d: signal/pedestal difference", i));

   sgn_diff->SetLineColor(kBlack);
   sgn_diff->SetFillColor(kBlue-10);
   sgn_diff->SetFillStyle(3001);
   //print_extra_info(sgn_diff);


	//superimpose histogramm
    sgn_diff->DrawCopy("same");
   

   

	/**
	* 7)
	* efficiency calculation and rescale
	**/
	efficiency = efficiency_calc(sgn_diff, sgn_tot_histogram);
	efficiency = efficiency	/ scale_factor;

   GGM_analysis_log( Form("Channel %i efficiency is %g\n\n", i, efficiency) );
   
   return efficiency;
}


/**
* Take a ntuple and return a histogram with range and binning calculated to remove outliers
*
* @param TNtuple
* @param Int_t index to take the correct branch (a comun in data file) from the ntuple
* @return new histogram width new binning, or NULL if something wrong
*/
TH1F *RemoveOutliers(TNtuple* ntupla, Int_t ch_number){

   TH1F* hist_temp = NULL;

   gDirectory->Delete("htemp1"); //clean memory
   ntupla->Draw( Form("channel_%d>>htemp1", ch_number),"","goff" );
   hist_temp = (TH1F*)gDirectory->FindObject("htemp1");
   
   	//return no histogram was drawn
   if( hist_temp == NULL){
      return NULL;
   } 
   

   Double_t *q = new Double_t[3];
   Double_t *p = new Double_t[3];
   q[0] = 0.; q[1] = 0.; q[2] = 0.;
   p[0] = 0.25; p[1] = 0.5; p[2] = 0.75;

   hist_temp->GetQuantiles(3,q,p);

   Double_t iqr = q[2] - q[0];
   Double_t up = TMath::Ceil( q[2] + 2*iqr );
   Double_t low = TMath::Floor( q[0] - 2*iqr );
   
   //TEST adaptable bin with in relation with the number of entries
   //Double_t bin_width = TMath::Ceil( (2*iqr)/TMath::Power( hist_temp->GetEntries(), 1/3.0 ) );

	Double_t bin_width = BIN_WIDTH;
   	up = RoundUp(up,bin_width);
	low = RoundUp(low,bin_width);
	

   Int_t bins = TMath::Nint( (up-low)/ bin_width );
	
	//return if no bins
   if( bins < 1){
      return NULL;
   }
	
	
   ntupla->Draw(Form("channel_%d>>htemp2(%d,%g,%g)", ch_number, bins, low, up) ,"", "goff");
   hist_temp = (TH1F*)gDirectory->FindObject("htemp2");

   
   return hist_temp;
}

/**
* Taken an histogram, find first bin with content <= 0
*
* @param TH1 histogram
* @return content of the first bin with content <= 0
**/
Double_t FindFirstZeroBeforeMaximun(TH1* hist){

   for (Int_t bin = hist->GetMaximumBin(); bin>=1; bin--) {
      if (hist->GetBinContent(bin) <= 0){
         return hist->GetBinCenter(bin);
      }

   }


   return hist->GetBinContent(0); //return the underflow
}

/**
* round up the first param to the next multiple of second parameters
*
* @param number to round up
* @param number to rounf with
* @return rounded number multiple of second param**/
Double_t RoundUp(Double_t i, Double_t n) {
    
    if(n > 0)
        return TMath::Ceil(i/n) * n;
    else if( n < 0)
        return TMath::Floor(i/n) * n;
    else
        return n;

}



/**
*  populate the NTuple with data from file in file_path param
* data are supposed to be 20 columns delimited by a whitespace 
*
* @param file_path, path to raw data file
* @param name, a name to assign to Ntuple
*
* @return pointer to TNtuple popolated with raw data from file path, NULL if error
**/
TNtuple* populate_tree(TString file_path, TString name){

   TNtuple *ntupla = (TNtuple*)gDirectory->FindObject(name); //se la NTupla è presente in memoria perchè già creata in precedenza, riutilizza la stessa senza processare di nuovo il file di input
   if( ntupla != NULL ){
	   gDirectory->Delete(name); //delete and recreate if already exists
   }
    
	ntupla = new TNtuple(name, "Data from ascii file", "event_id/i:column_2/D:timestamp:column_4:channel_1/D:channel_2:channel_3:channel_4:channel_5:channel_6:channel_7:channel_8:channel_9:channel_10:channel_11:channel_12:channel_13:channel_14:channel_15:channel_16");
      
	if( !(ntupla->ReadFile(file_path)) ){
		return NULL;
   }
   
	return ntupla;
}

/**
*
* stampa informazioni sull'istogramma passato come argomento
*
* @param hist puntatore a TH1
**/
void print_extra_info(TH1 *hist){

   Double_t ch_mean = 0;
   Double_t ch_rms = 0;
   Double_t lower_limit = 0;
   Double_t upper_limit = 0;

   TAxis *hist_xaxis = hist->GetXaxis(); //prendo asse X istogramma

   ch_mean = hist->GetMean(); //media istogramma
   ch_rms = hist->GetRMS(); //root mean square

   lower_limit = hist_xaxis->GetXmin();
   upper_limit = hist_xaxis->GetXmax();


   cout << hist->GetTitle() << endl;
   cout << "  number of entries: " << hist->GetEntries() << endl;
   cout << "  Mean. " << ch_mean << endl;
   cout << "  RMS. " << ch_rms << endl;
   cout << "  lower_limit. " << lower_limit << endl;
   cout << "  upper_limit. " << upper_limit << endl;
   cout << "  number of bins: " << hist_xaxis->GetNbins() << endl;
   cout << "  bin width: " << hist_xaxis->GetBinWidth(hist->GetMaximumBin()) <<  endl;
   cout << "  max: " << hist_xaxis->GetBinCenter(hist->GetMaximumBin()) << endl << endl;

    //Controllo esistenza di un canvas attivo
   if( gPad != NULL ){
      gPad->Update(); //necessario per forza la creazione del box stats
      TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
      if( ps != NULL ){ //verifico se il box stats esiste e ci aggiungo due informazioni in più
            ps->SetName("mystats");
            //TList *list = ps->GetListOfLines();

            ps->AddText(Form("bin width = %g",(Double_t)hist_xaxis->GetBinWidth(hist->GetMaximumBin())));
            ps->AddText(Form("number of bins = %g", (Double_t)hist_xaxis->GetNbins()));

            // the following line is needed to avoid that the automatic redrawing of stats
            hist->SetStats(0);
      }

      
   gPad->Modified();
   }

}


/**
*
* @param1 TH1*, Real signal histogram (results from the difference between total signal and pedestal)
* @param2 TH1*, Total signal histogtam
*
* @return chamber efficiency
*
**/
Double_t efficiency_calc(TH1 *sgn_diff, TH1 *sgn_tot){

   
   // integral from the bin where x=0 to the last bin+1 to include overflows
   Int_t fromBin = sgn_diff->GetXaxis()->FindFixBin(0);
   Int_t toBin = sgn_diff->GetXaxis()->GetLast() + 1;
   Double_t sgn_diff_area = sgn_diff->Integral(fromBin, toBin);
 
	
	GGM_analysis_log( Form("Real signal area: %g\n", sgn_diff_area) );
	
	
	// full range integral equals number of entries
	// half because signal histogram as double entries
	Double_t sgn_tot_area = sgn_tot->GetEntries()/2; 
	
	GGM_analysis_log( Form("Total signal area: %g\n\n", sgn_tot_area/sgn_diff_area) );

	
   return sgn_diff_area/sgn_tot_area; // ratio area is the efficiency
}


/**
*
* crea il file .DST contenente
* 1. timestamp
* 2. numero canale
* 3 a 6. uno zero
* 7. efficenza
*
* @param index, indice canale
* @param efficiency,  l'efficenza del canale calcolata
* @param filename, nome del file in cui salvare i dati
* @return void
*
**/
void update_dst_file(int index, double efficiency, TString filename){
	
	TString output_filename;
	
	if( filename.IsNull() )
		output_filename = "temp.dst";
	else
		output_filename = Form("%s.dst", filename.Data());
	
	
   ofstream myFile(output_filename.Data(), ofstream::app); //apre il file in append per aggiungere righe al file

   TTimeStamp *timestamp = new TTimeStamp();
   
   //myFile << Form("%u %d 0 0 0 0 %g\n", timestamp->GetSec(), index, efficiency);
   /*** TEST ***/
	double hv_value = read_hv(Form("./DSToriginali/%s_l_3.dst",filename.Data()), index-1); //leggo dal file originale il valore di HV, il contatore dei canali inizia da zero
	
   myFile << Form("%u %d 0 0 0 0 %g %g\n", timestamp->GetSec(), index, efficiency, hv_value);
   /************/
   //myFile << Form("%u %d 0 0 0 0 %g\n", timestamp->GetSec(), index, efficiency);

   myFile.close();

	return;
}

/**
*
* fa un loop tra tutti i bin e imposta a zero quelli con valori negativi
*
* @param puntatore a TH1
* @return void
*  
*/
void set_negative_to_zero(TH1* hist){
   Int_t nbinsx = 0;
   Int_t bin_content=0;

   nbinsx = hist->GetNbinsX();

   for(Int_t i=1; i <= nbinsx; i++){
      if( hist->GetBinContent(i) < 0){
         hist->SetBinContent(i,0);
      }
   }
}

/**
*
* fa un loop tra tutti i bin e imposta a zero quelli sulle ascisse negative
*
* @param puntatore a TH1
* @return void
*  
*/
void set_negative_bin_to_zero(TH1* hist){

   Int_t nbinsx = 0;
   Int_t bin_content=0;

   TAxis* asseX = hist->GetXaxis();
   nbinsx = asseX->GetNbins();

   for(Int_t i=1; i <= nbinsx; i++){
	   if( asseX->GetBinLowEdge(i) < 0 ){
		   hist->SetBinContent(i,0);
	   }	
   }
}



double read_hv(TString filename, int k){
	
	Int_t nlines = 0;
	double a,b,c,d,e,f,g,h,i,l; //colonna b = indice canale, colonna h = HV
	
	ifstream in;
	in.open(filename.Data());
	
	 while (1) {
		  in >> a >> b >> c >> d >> e >> f >> g >> h >> i >> l;
		  if (!in.good())
			  break;
		  if( (double)k == b ){
			  printf("canale: %d\nHV letta: %g\n\n\n",b+1,h);
			  return h;
		  }
		  
		  nlines++;
   }
   
   //printf("canale: %d non presente nel dst\n\n\n",k+1);
   return 0;
}

/**
*
* take an histogram and set label on his axis
*
* @param TH1*
* @param TString label for abscissa
* @param TString label for ordinate
*
**/
void setAxisTitle(TH1* hist, TString Xlabel, TString Ylabel){
	
	TAxis* asseX = hist->GetXaxis();
	asseX->SetTitle(Xlabel);
	asseX->SetTitleOffset(1.4);
	asseX->SetLabelSize(0.030);
	
	TAxis* asseY = hist->GetYaxis();
	asseY->SetTitle(Ylabel);
	asseY->SetTitleOffset(1.6);
	asseY->SetLabelSize(0.030);
	
	return;
}