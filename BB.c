void BB();
void exit_here();
void particle();
void atomic_number();
double bethe_bloch();
void draw_graph();
void plot_dissipazione();



TCanvas *main_window;
TButton *exit_button; 
TButton *alpha_button;
TButton *proton_button;
TButton *muon_button;
TButton *graph_button;
TButton *gold_button;
TButton *al_button;
TButton *alluminio_button;
TButton *fe_button;
TButton *dissipazione_button;
TButton *clear_button;
TPad *pad1;
TPad *pad2;

//Inizializzo le variabili
double m; //massa del proiettile
int Z; //numero atomico del bersaglio
float A; //massa atomica relativa
float z_particle; //numero atomico proiettile
double mrel; //

//Inizializzo e dichiaro le costanti globali
double c = 2.9979*pow(10,8); //velocità luce in m/s
double me = 0.511*pow(10,6); //massa a riposo elettrone in eV/c^2
int e = 1; //carica elettrone in eV
double K = 0.307075*pow(10,6); //costante K
long int incrementi = 10000; //step di integrazione e di calcolo dell'equazione
float betagamma_max = 10000; 
float betagamma_min = 0.1;
float TEXT_SIZE=0.30;
float MAX_SPESSORE = 1000;

//valori fissi per l'interfaccia grafica
Double_t w = 1250;
Double_t h = 800;

void BB() 
{
	TCanvas * main_window = new TCanvas("main_window","Bethe-Bloch");
	TPad* pad1 = new TPad("pad1", "Graph", 0, 0, 0.9, 1);
	TPad* pad2 = new TPad("pad2", "Buttons", 0.90, 0, 1, 1, 38);
	pad1->Draw();
	pad2->Draw();
	main_window->SetCanvasSize(w, h);
	main_window->SetWindowSize(w, h);
	
	//Uso pad2 per organizzare i miei bottoni
	pad2->cd();
	
	//Creo dei bottoni che aiutino l'utente nella scelta della configurazione iniziale


	proton_button = new TButton("Protone", "particle(0)", 0, 0.95, 1, 0.90);
	proton_button->SetTextSize(TEXT_SIZE);
	proton_button->Draw();

	alpha_button = new TButton("Particella Alfa", "particle(1)", 0,0.9,1,0.85);	
	alpha_button->SetTextSize(TEXT_SIZE);
	alpha_button->Draw();
	
	muon_button = new TButton("Muone", "particle(2)", 0,0.85,1,0.8);	
	muon_button->SetTextSize(TEXT_SIZE);
	muon_button->Draw();
	
	//Bottone "ferro"
	fe_button = new TButton("Ferro", "atomic_number(0)", 0, 0.75, 1, 0.7);
	fe_button->SetTextSize(TEXT_SIZE);
	fe_button->Draw();

	//Bottone "oro"
	gold_button = new TButton("Oro", "atomic_number(1)", 0, 0.7, 1, 0.65);
	gold_button->SetTextSize(TEXT_SIZE);
	gold_button->Draw();

	//Bottone "argento"
	al_button = new TButton("Argento", "atomic_number(2)", 0, 0.65, 1, 0.6);
	al_button->SetTextSize(TEXT_SIZE);
	al_button->Draw();
	
	//bottone "alluminio"
	alluminio_button = new TButton("Alluminio", "atomic_number(3)", 0, 0.6, 1, 0.55);
	alluminio_button->SetTextSize(TEXT_SIZE);
	alluminio_button->Draw();


	graph_button = new TButton("Grafico Bethe-Bloch", "draw_graph()", 0, 0.4,1, 0.35);
	graph_button->SetTextSize(TEXT_SIZE);
	graph_button->Draw();

	dissipazione_button = new TButton("Grafico dissipazione","plot_dissipazione()", 0, 0.35,1,0.3);
	dissipazione_button->SetTextSize(TEXT_SIZE);
	dissipazione_button->Draw();
	
	clear_button = new TButton("Cancella grafico", "", 0, 0.3,1,0.25);
	clear_button->SetTextSize(TEXT_SIZE);
	clear_button->Draw();
	
	exit_button = new TButton("Esci", "exit_here()", 0, 0.1,1, 0.15);
	exit_button->SetTextSize(TEXT_SIZE);
	exit_button->Draw();
	pad1->cd();
}
//con questo codice chiudo il software
void exit_here() 
{
        cout << "exiting root by running\ngROOT->ProcessLine(\".q\")\n";
        gROOT->ProcessLine(".q");
}

//Qui puoi selezionare il tuo proiettile e caricarne i valori
void particle(int k)
{
	if (k == 0){
		m = 9.38257*pow(10,8);
		z_particle = 1;
		cout << "Protone selezionato." << endl;
		} 
	else if (k == 1) { 
		m = 3.7273*pow(10,9);
		z_particle = 2;
		cout << "Particella alfa selezionata." << endl; }
	else if (k == 2) {
		m = 1.0565*pow(10,8);
		z_particle = 1; 
		cout << "Muone selezionato." << endl;}
}
		
//Qui puoi selezionare il bersaglio
void atomic_number(int l)
{
	if (l == 0){
		Z = 26;
		A = 53.9396;
		cout << "Ferro selezionato." << endl;
		} 
	if (l == 1) { 
		Z = 79;
		A = 196.9666;
		cout << "Oro selezionato." << endl; 
		}
	if (l == 2) { 
		Z = 47;
		A = 106.9051;
		cout << "Argento selezionato." << endl; 
		}
	if (l == 3){
		Z = 13;
		A = 26.9815;
		cout << "Alluminio selezionato." << endl;
		
	}
}

	//Bethe-bloch function
double bethe_bloch(double gam, double beta)
{
	mrel = me/m;
	double Tmax = 2*me*pow(gam*beta,2)/(1+(2*gam*mrel)+pow(mrel,2));
	//determino I
	double I = Z*10;
	//applico Bethe-bloch 
	double func = K*pow(z_particle,2)*(Z/A)*(1/pow(beta,2))*(0.5*log(2*me*pow(beta*gam,2)*Tmax/pow(I,2))-pow(beta,2));
	//restituisco il valore
	return func;
}

//Plotto i valori
void draw_graph()
{
	double beta_min = betagamma_min/sqrt(1+pow(betagamma_min,2));
	double beta_max = betagamma_max/sqrt(1+pow(betagamma_max,2));
	TGraph* graph = new TGraph(); //creo un grafico
	graph->SetTitle("Formula di Bethe-Bloch; #gamma*#beta; - #frac{dE}{dx} [MeV]"); //imposto il titolo del grafico e il nome degli assi
   	graph->GetXaxis()->CenterTitle(true);	
   	graph->GetYaxis()->CenterTitle(true);	
		//Assegno diversi colori in base al materiale del bersaglio
		if (Z == 26){
		graph->SetLineColor(28);
		}
		if (Z == 47){
		graph->SetLineColor(16);
		}
		if (Z == 79){
		graph->SetLineColor(5);
		}
		if (Z == 13){
		graph->SetLineColor(12);
		}
	if (Z == 0 || m == 0){
		cout << "Seleziona una configurazione iniziale." << endl;
		}
	else {
		double x,y;
		double delta_beta = (beta_max-beta_min)/incrementi;
		long int i=0;
		for (double beta = beta_min; beta < beta_max; beta = beta+delta_beta) {
			//gamma
			double gam = 1/sqrt(1-beta*beta);
			y=bethe_bloch(gam, beta)*pow(10,-6);
			x=beta*gam;
			graph->SetPoint(i++,x,y);
		}
	//applico alcune correzioni stilistiche e restituisco il grafico
	graph->SetLineWidth(4);
	gPad->SetLogx();
	gPad->SetLogy();
	graph->Draw();
	}
	
	
}
void plot_dissipazione()
	{
	if (Z == 0 || m == 0){
	cout << "Seleziona una configurazione iniziale." << endl;
	}
	else {
		
		TGraph* graph = new TGraph(); //creo un grafico
		graph->SetTitle("Dissipazione Energia; Lunghezza percorso [cm]; - Stopping Power #frac{[MeV]}{[cm]}");
		graph->GetXaxis()->CenterTitle(true);	
   		graph->GetYaxis()->CenterTitle(true);
		long int i=0;
		double x,inc, beta, gam, y,de;
		x=0;
		beta=0.6;
		gam = 1/sqrt(1-beta*beta);
		y=m*gam;
		de=bethe_bloch(gam,beta);
		inc = 0.005;
		for(;y>0;)
		{
		de=-bethe_bloch(gam,beta)*inc;
		graph->SetPoint(i++,x,-de*pow(10,-6));
		x = x+inc;
		y = y+de;
		beta = sqrt(1-((m*m)/(y*y)));
		gam = 1/sqrt(1-(beta*beta));
		}
		graph->SetLineWidth(4);
		gPad->SetLogy();
		//gPad->SetLogx();
		graph->Draw();
	}
	}
