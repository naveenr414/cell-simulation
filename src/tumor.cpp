/* 

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
#include <stdio.h>
#include <malloc.h>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
//#include <unistd.h>
#include <math.h>
#include "dish.h"
#include "random.h"
#include "cell.h"
#include "info.h"
#include "parameter.h"
#include "sqr.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif


using namespace std;

INIT {

  try {
    
    // Define initial distribution of cells
    //CPM->GrowInCells(1, 10, 1); // initial tumor cell
    //CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells, (par.sizex/2-2), par.sizey-2, par.sizex/2 ,0);
    CPM->ConstructInitCells(*this);
    
    // Set the first cell to type 2
    /* getCell(1).setTau(2);
    getCell(1).SetTargetLength(10);
    getCell(1).SetTargetArea(50);*/
    
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }

}

TIMESTEP { 
 
  try {

    static int i=0;
  
    static Dish *dish=new Dish();
    static Info *info=new Info(*dish, *this);

    // slowly increase target length during the first time steps
    // to prevent cells from breaking apart
    // static double targetlength=par.target_length;


     for (int r=0;r<par.pde_its;r++) {
	
	dish->PDEfield->Secrete(dish->CPM);
	dish->PDEfield->Diffuse(1);

      }
      dish->CPM->AmoebaeMove(dish->PDEfield);
     
    if (i>=par.relaxation) {
      
      if (!(i%5)) {
	dish->MeasureChemConcentrations();
	dish->CPM->GrowAndDivideCells(2);
      }
      
    }
    // Plot
    if (par.graphics && !(i%par.storage_stride)) {
      
      int tx,ty;
      
      BeginScene();
            
     
      dish->PDEfield->Plot(this,0);
      
      // You need to call "ClearImage" if no PDE field is plotted,
      // because the CPM medium is considered transparant
      //ClearImage();
      dish->Plot(this);
 
      
      if (i>=par.relaxation)
	 dish->PDEfield->ContourPlot(this,0,7);
      
      char title[400];
      snprintf(title,399,"CellularPotts: %.2f hr",dish->PDEfield->TheTime()/3600);
      //snprintf(title,399,"CellularPotts: %d MCS",i);
      //ChangeTitle(title);
      EndScene();
      info->Menu();
          
    }
    
    if (par.store && !(i%par.storage_stride)) {

      char fname[200];
      sprintf(fname,"%s/extend%05d.png",par.datadir,i);
    
      BeginScene();

      dish->PDEfield->Plot(this,0);
      //ClearImage();
      dish->Plot(this);
      if (i>=par.relaxation)
	dish->PDEfield->ContourPlot(this,0,7);
   
      EndScene();
      
      Write(fname);
        
    }

    i++;
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }
}

void PDE::Secrete(CellularPotts *cpm) {

  const double dt=par.dt;
  int ctype;
  
  for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++) {
      // inside cells
      if (ctype=cpm->Sigma(x,y)) {
	
	sigma[1][x][y]-=par.decay_rate[1]*sigma[1][x][y]*dt;
	if (cpm->getCell(ctype).getTau()!=2) { 
	  sigma[0][x][y]+=dt*(par.secr_rate[0]);
	}
      } else {
      // outside cells (source term)
	sigma[1][x][y]=1.;
	sigma[0][x][y]-=dt*(par.decay_rate[0]*sigma[0][x][y]);
      }
    }
}

int PDE::MapColour(double val) {
  
  return (((int)((val/((val)+1.))*100))%100)+155;
}

int main(int argc, char *argv[]) {
  
	
  try {

#ifdef QTGRAPHICS
    QApplication a(argc, argv);
#endif
    // Read parameters
    par.Read(argv[1]);
    
    Seed(par.rseed);
    
    //QMainWindow mainwindow w;
#ifdef QTGRAPHICS
    QtGraphics g(par.sizex*2,par.sizey*2);
    a.setMainWidget( &g );
    a.connect(&g, SIGNAL(SimulationDone(void)), SLOT(quit(void)) );

    if (par.graphics)
      g.show();
    
    a.exec();
#else
    X11Graphics g(par.sizex*2,par.sizey*2);
    int t;

    for (t=0;t<par.mcs;t++) {

      g.TimeStep();
    
    }
#endif
    
  } catch(const char* error) {
    std::cerr << error << "\n";
    exit(1);
  }
  catch(...) {
    std::cerr << "An unknown exception was caught\n";
  }
  return 0;
}
