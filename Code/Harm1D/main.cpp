/**                                                                                        
 * GENERAL REMARKS                                                                         
 *                                                                                         
 *  This code is freely available under the following conditions:                          
 *                                                                                         
 *  1) The code is to be used only for non-commercial purposes.                            
 *  2) No changes and modifications to the code without prior permission of the developer. 
 *  3) No forwarding the code to a third party without prior permission of the developer.  
 *                                                                                         
 *  			MTCalc_with_DFP_COCR                                               
 *   Program entry for 1D MT problem    
 *                                                                                         
 *  Written by Ph.D. Petr A. Domnikov                                                      
 *  Novosibirsk State Technical University,                                                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                      
 *  p_domnikov@mail.ru                                                                     
 *  Version 1.2 December 5, 2020                                                          
*/                                                                                         


#include "stdafx.h"
#include "divgrad_1d.h"

ofstream logfile;
//------------------------------------------------------------
int main()
{
	logfile.open("LogMTZ1D");

	Divgrad_1d divgrad1d;

	divgrad1d.Solve_1d_Problem_for_3d_task();

	logfile.close();
	logfile.clear();

	return 0;
}