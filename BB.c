//Funzioni che ci serviranno man mano
void BB();
void exit_here();
void particle();
void atomic_number();
double bethe_bloch();
void draw_graph();
void plot_dissipazione();

//Puntatori relativi ad elementi grafici
TCanvas *main_window;
TButton *exit_button; 
TButton *alpha_button;
TButton *proton_button;
TButton *muon_button;
TButton *graph_button;
TButton *au_button;
TButton *ag_button;
TButton *al_button;
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
double mrel; //massa relativa del proiettile rispetto a quella dell'elettrone

//Inizializzo e dichiaro le costanti globali
double c = 2.9979*pow(10,8); //velocità luce in m/s
double me = 0.511*pow(10,6); //massa a riposo elettrone in eV/c^2
double K = 0.307075*pow(10,6); //costante K

//Inizializzo e dichiaro le costanti selezionabili dall'utente
long int INCREMENTI = 10000; //step di calcolo dell'equazione
float BETAGAMMA_MAX = 10000; //valore massimo di beta*gamma che compare nel grafico di Bethe-Bloch
float BETAGAMMA_MIN = 0.1; //valore minimo
float TEXT_SIZE=0.30; //grandezza del font dei bottoni
float INC = 0.001; //spessore infinitesimo in cm;

void BB()  //Interfaccia utente
{
	TCanvas * main_window = new TCanvas("main_window","Bethe-Bloch");
	TPad* pad1 = new TPad("pad1", "Graph", 0, 0, 0.9, 1);
	TPad* pad2 = new TPad("pad2", "Buttons", 0.90, 0, 1, 1, 38);
	pad1->Draw();
	pad2->Draw();

	//In pad2 andranno i bottoni
	pad2->cd();
	
	//Creo dei bottoni che aiutino l'utente nella scelta del proiettile...
	proton_button = new TButton("Protone", "particle(0)", 0, 0.95, 1, 0.90);
	proton_button->SetTextSize(TEXT_SIZE);
	proton_button->Draw();

	alpha_button = new TButton("Particella Alfa", "particle(1)", 0,0.9,1,0.85);	
	alpha_button->SetTextSize(TEXT_SIZE);
	alpha_button->Draw();
	
	muon_button = new TButton("Muone", "particle(2)", 0,0.85,1,0.8);	
	muon_button->SetTextSize(TEXT_SIZE);
	muon_button->Draw();
	
	//... e del bersaglio
	fe_button = new TButton("Ferro", "atomic_number(0)", 0, 0.75, 1, 0.7);
	fe_button->SetTextSize(TEXT_SIZE);
	fe_button->Draw();

	au_button = new TButton("Oro", "atomic_number(1)", 0, 0.7, 1, 0.65);
	au_button->SetTextSize(TEXT_SIZE);
	au_button->Draw();

	ag_button = new TButton("Argento", "atomic_number(2)", 0, 0.65, 1, 0.6);
	ag_button->SetTextSize(TEXT_SIZE);
	ag_button->Draw();
	
	al_button = new TButton("Alluminio", "atomic_number(3)", 0, 0.6, 1, 0.55);
	al_button->SetTextSize(TEXT_SIZE);
	al_button->Draw();

	//Bottoni che richiamano funzioni ed eseguono comandi:
	graph_button = new TButton("Grafico Bethe-Bloch", "draw_graph()", 0, 0.4,1, 0.35);
	graph_button->SetTextSize(TEXT_SIZE);
	graph_button->Draw();

	dissipazione_button = new TButton("Grafico dissipazione","plot_dissipazione()", 0, 0.35,1,0.3);
	dissipazione_button->SetTextSize(TEXT_SIZE);
	dissipazione_button->Draw();
	
	clear_button = new TButton("Cancella grafico", ")", 0, 0.3,1,0.25); //DA FARE
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

void particle(int k) //Qui puoi selezionare il tuo proiettile e caricarne i valori
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
		

void atomic_number(int l) //Qui puoi selezionare il bersaglio
{
	if (l == 0) {
		Z = 26;
		A = 53.9396;
		cout << "Ferro selezionato." << endl;} 
	if (l == 1) { 
		Z = 79;
		A = 196.9665;
		cout << "Oro selezionato." << endl; }
	if (l == 2) { 
		Z = 47;
		A = 106.9051;
		cout << "Argento selezionato." << endl; }
	if (l == 3){
		Z = 13;
		A = 26.9815;
		cout << "Alluminio selezionato." << endl;}
}

double bethe_bloch(double beta) //Uno strumento che utilizzeremo molto spesso
{
	mrel = me/m;
	double gam = 1/sqrt(1-pow(beta,2));
	double Tmax = 2*me*pow(gam*beta,2)/(1+(2*gam*mrel)+pow(mrel,2));
	//determino I
	double I = Z*10;
	//applico Bethe-bloch 
	double func = K*pow(z_particle,2)*(Z/A)*pow(beta,-2)*(0.5*log(2*me*Tmax*pow(beta*gam/I,2))-pow(beta,2));
	//restituisco il valore
	return func;
}

void draw_graph() //Calcolo e plotto i valori di Bethe_Bloch
{
	double beta_min = BETAGAMMA_MIN/sqrt(1+pow(BETAGAMMA_MIN,2));
	double beta_max = BETAGAMMA_MAX/sqrt(1+pow(BETAGAMMA_MAX,2));
	TGraph* graph = new TGraph(); //creo un grafico
	graph->SetTitle("Formula di Bethe-Bloch; #gamma*#beta; - #frac{dE}{dx} [MeV]"); //imposto il titolo del grafico e il nome degli assi
   	graph->GetXaxis()->CenterTitle(true);	
   	graph->GetYaxis()->CenterTitle(true);	
	
	//Assegno diversi colori in base al materiale del bersaglio
	if (Z == 26){graph->SetLineColor(28);}
	if (Z == 47){graph->SetLineColor(16);}
	if (Z == 79){graph->SetLineColor(5);}
	if (Z == 13){graph->SetLineColor(12);}
	if (Z == 0 || m == 0){cout << "Seleziona una configurazione iniziale." << endl;}
	else 
	{
		double x,y;
		double delta_beta = (beta_max-beta_min)/INCREMENTI;
		long int i=0;
		for (double beta = beta_min; beta < beta_max; beta = beta+delta_beta)
		{
			double gam = 1/sqrt(1-beta*beta);
			y=bethe_bloch(beta)*pow(10,-6); //con un fattore di scala 10^-6 per mostrare i valori in MeV
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
		
void plot_dissipazione() //mostro il grafico della Stopping Power fino a quando la particella si arresta
{
	if (Z == 0 || m == 0)
	{
		cout << "Seleziona una configurazione iniziale." << endl;
	}
	else
	{
		TGraph* graph = new TGraph(); //creo un grafico
		graph->SetTitle("Dissipazione Energia; Spessore percorso dalla particella [cm]; Stopping Power #frac{[MeV]}{[cm]}");
		graph->GetXaxis()->CenterTitle(true);	
   		graph->GetYaxis()->CenterTitle(true);
		double x, beta, gam, y,de;
		long int i=0;
		x=0; //distanza percorsa nel bersaglio
		beta=0.2; //fornisco il beta della particella
		gam = 1/sqrt(1-pow(beta,2)); //calcolo il corrispondente gamma
		y=m*gam; //energia cinetica relativistica della particella prima che impatti contro il materiale
		for(;y>=0;i++) //eseguo il ciclo finchè l'energia cinetica è 0
		{
			de=-bethe_bloch(beta)*INC; //perdita di energia infinitesima
			graph->SetPoint(i,x,-de*pow(10,-6)); //faccio il plot "spazio percorso vs energia dissipata"
			x = x+INC; //incremento la distanza percorsa nel bersaglio di un determinato spessore infinitesimo
			y = y+de; //calcolo la perdita di energia cinetica (col segno più perchè de è negativa)
			beta = sqrt(1-(pow(m/y,2))); //calcolo il nuovo beta corrispondente alla nuova energia
			gam = 1/sqrt(1-(pow(beta,2))); //calcolo il nuovo gamma corrispondente al nuovo beta
		}
		graph->SetLineWidth(3);
		gPad->SetLogy();
		//gPad->SetLogx();
		graph->Draw();
	}
}	
