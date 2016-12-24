//////////////////////////
// parser.hpp
//////////////////////////
// 
// Parser for settings file
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: November 2016
//
//////////////////////////

#ifndef PARSER_HEADER
#define PARSER_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "metadata.hpp"

using namespace std;

struct parameter
{
	char name[PARAM_MAX_LENGTH];
	char value[PARAM_MAX_LENGTH];
	bool used;
};


int sort_descending(const void * z1, const void * z2)
{
	return (* (double *) z1 > * (double *) z2) ? -1 : ((* (double *) z1 < * (double *) z2) ? 1 : 0);
}


//////////////////////////
// readline
//////////////////////////
// Description:
//   reads a line of characters and checks if it declares a parameter; if yes,
//   i.e. the line has the format "<parameter name> = <parameter value>" (with
//   an optional comment added, preceded by a hash-symbol, '#'), the parameter
//   name and value are copied to the corresponding arrays, and 'true' is returned.
//   If the format is not recognized (or the line is commented using the hash-symbol)
//   'false' is returned instead.
// 
// Arguments:
//   line       string containing the line to be read
//   pname      will contain the name of the declared parameter (if found)
//   pvalue     will contain the value of the declared parameter (if found)
//
// Returns:
//   'true' if a parameter is declared in the line, 'false' otherwise.
// 
//////////////////////////

bool readline(char * line, char * pname, char * pvalue)
{
	char * pequal;
	char * phash;
	char * l;
	char * r;
	
	pequal = strchr(line, '=');
	
	if (pequal == NULL || pequal == line) return false;
	
	phash = strchr(line, '#');
	
	if (phash != NULL && phash < pequal) return false;
	
	l = line;
	while (*l == ' ' || *l == '\t') l++;
	
	r = pequal-1;
	while ((*r == ' ' || *r == '\t') && r > line) r--;
	
	if (r < l) return false;
	
	if (r-l+1 >= PARAM_MAX_LENGTH) return false;
	
	strncpy(pname, l, r-l+1);
	pname[r-l+1] = '\0';
	
	l = pequal+1;
	while (*l == ' ' || *l == '\t') l++;
	
	if (phash == NULL)
		r = line+strlen(line)-1;
	else
		r = phash-1;
	  
	while (*r == ' ' || *r == '\t' || *r == '\n' || *r == '\r') r--;
	
	if (r < l) return false;
	
	if (r-l+1 >= PARAM_MAX_LENGTH) return false;
	
	strncpy(pvalue, l, r-l+1);
	pvalue[r-l+1] = '\0';
	
	return true;
}


//////////////////////////
// loadParameterFile
//////////////////////////
// Description:
//   loads a parameter file and creates an array of parameters declared therein
// 
// Arguments:
//   filename   string containing the path to the parameter file
//   params     will contain the array of parameters (memory will be allocated)
//
// Returns:
//   number of parameters defined in the parameter file (= length of parameter array)
// 
//////////////////////////

int loadParameterFile(const char * filename, parameter * & params)
{
	int numparam = 0;
	int i = 0;
	
	if (parallel.grid_rank()[0] == 0) // read file
	{
		FILE * paramfile;
		char line[PARAM_MAX_LINESIZE];
		char pname[PARAM_MAX_LENGTH];
		char pvalue[PARAM_MAX_LENGTH];
		
		paramfile = fopen(filename, "r");
		
		if (paramfile == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Unable to open parameter file " << filename << "." << endl;
			parallel.abortForce();
		}
		
		while (!feof(paramfile) && !ferror(paramfile))
		{
			fgets(line, PARAM_MAX_LINESIZE, paramfile);
			
			if (readline(line, pname, pvalue) == true) numparam++;
		}
		
		if (numparam == 0)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! No valid data found in file " << filename << "." << endl;
			fclose(paramfile);
			parallel.abortForce();
		}
		
		params = (parameter *) malloc(sizeof(parameter) * numparam);
		
		if (params == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Memory error." << endl;
			fclose(paramfile);
			parallel.abortForce();
		}
		
		rewind(paramfile);
		
		while (!feof(paramfile) && !ferror(paramfile) && i < numparam)
		{
			fgets(line, PARAM_MAX_LINESIZE, paramfile);
			
			if (readline(line, params[i].name, params[i].value) == true)
			{
				params[i].used = false;
				i++;
			}
		}
		
		fclose(paramfile);
		
		if (i < numparam)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! File may have changed or file pointer corrupted." << endl;
			free(params);
			parallel.abortForce();
		}
		
		parallel.broadcast_dim0<int>(numparam, 0);
	}
	else
	{
		parallel.broadcast_dim0<int>(numparam, 0);
		
		params = (parameter *) malloc(sizeof(parameter) * numparam);
		
		if (params == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Memory error." << endl;
			parallel.abortForce();
		}
	}
	
	parallel.broadcast_dim0<parameter>(params, numparam, 0);
	
	return numparam;
}


//////////////////////////
// saveParameterFile
//////////////////////////
// Description:
//   saves a parameter file
// 
// Arguments:
//   filename   string containing the path to the parameter file
//   params     array of parameters
//   numparam   length of parameter array
//   used_only  if 'true', only the used parameters will be written (default)
//
// Returns:
// 
//////////////////////////

void saveParameterFile(const char * filename, parameter * params, const int numparam, bool used_only = true)
{
	if (parallel.isRoot())
	{
		FILE * paramfile;
		
		paramfile = fopen(filename, "w");
		
		if (paramfile == NULL)
		{
			cout << " error in saveParameterFile! Unable to open file " << filename << "." << endl;
		}
		else
		{
			for (int i = 0; i < numparam; i++)
			{
				if (!used_only || params[i].used)
					fprintf(paramfile, "%s = %s\n", params[i].name, params[i].value);
			}
			
			fclose(paramfile);
		}
	}
}


//////////////////////////
// parseParameter (int)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as integer
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, int & pvalue)
{   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			if (sscanf(params[i].value, "%d", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}
	
	return false;
}


//////////////////////////
// parseParameter (long)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as integer
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, long & pvalue)
{   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			if (sscanf(params[i].value, "%ld", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}
	
	return false;
}


//////////////////////////
// parseParameter (double)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as double
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to double which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, double & pvalue)
{   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			if (sscanf(params[i].value, "%lf", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}
	
	return false;
}


//////////////////////////
// parseParameter (char *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and retrieves its value as string
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     character string which will contain a copy of the parameter value (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, char * pvalue)
{   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			strcpy(pvalue, params[i].value);
			params[i].used = true;
			return true;
		}
	}
	
	return false;
}


//////////////////////////
// parseParameter (double *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated double values
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of double values which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, double * pvalue, int & nmax)
{
	char * start;
	char * comma;
	char item[PARAM_MAX_LENGTH];
	int n = 0;
	   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if (nmax > 1)
			{
				while ((comma = strchr(start, ',')) != NULL)
				{
					strncpy(item, start, comma-start);
					item[comma-start] = '\0';
					if (sscanf(item, " %lf ", pvalue+n) != 1)
					{
						nmax = n;
						return false;
					}
					start = comma+1;
					if (++n > nmax-2)
						break;
				}
			}   
			if (sscanf(start, " %lf ", pvalue+n) != 1)
			{
				nmax = n;
				return false;
			}
			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}
	
	nmax = 0;
	return false;
}


//////////////////////////
// parseParameter (int *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated integer values
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of integer values which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, int * pvalue, int & nmax)
{
	char * start;
	char * comma;
	char item[PARAM_MAX_LENGTH];
	int n = 0;
	   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if (nmax > 1)
			{
				while ((comma = strchr(start, ',')) != NULL)
				{
					strncpy(item, start, comma-start);
					item[comma-start] = '\0';
					if (sscanf(item, " %d ", pvalue+n) != 1)
					{
						nmax = n;
						return false;
					}
					start = comma+1;
					if (++n > nmax-2)
						break;
				}
			}   
			if (sscanf(start, " %d ", pvalue+n) != 1)
			{
				nmax = n;
				return false;
			}
			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}
	
	nmax = 0;
	return false;
}


//////////////////////////
// parseParameter (char **)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated strings
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of character strings which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, char ** pvalue, int & nmax)
{
	char * start;
	char * comma;
	char * l;
	char * r;
	char item[PARAM_MAX_LENGTH];
	int n = 0;
	   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if (nmax > 1)
			{
				while ((comma = strchr(start, ',')) != NULL)
				{
					l = start;
					while (*l == ' ' || *l == '\t') l++;
					r = comma-1;
					while ((*r == ' ' || *r == '\t') && r > start) r--;
					
					if (r < l)
					{
						nmax = n;
						return false;
					}
					
					strncpy(pvalue[n], l, r-l+1);
					pvalue[n][r-l+1] = '\0';
					
					start = comma+1;
					if (++n > nmax-2)
						break;
				}
			}
			l = start;
			while (*l == ' ' || *l == '\t') l++;
			r = l;
			while (*r != ' ' && *r != '\t' && *r != '\0') r++;
			r--;
			
			if (r < l)
			{
				nmax = n;
				return false;
			}
			
			strncpy(pvalue[n], l, r-l+1);
			pvalue[n][r-l+1] = '\0';
			
			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}
	
	nmax = 0;
	return false;
}


//////////////////////////
// parseFieldSpecifiers
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated field specifiers
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     integer which will contain the binary-encoded list of parsed specifiers (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
// 
//////////////////////////

bool parseFieldSpecifiers(parameter * & params, const int numparam, const char * pname, int & pvalue)
{
	char * start;
	char * comma;
	int pos;
	char item[PARAM_MAX_LENGTH];
		   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			pvalue = 0;
			start = params[i].value;
			while ((comma = strchr(start, ',')) != NULL)
			{
				strncpy(item, start, comma-start);
				for (pos = comma-start; pos > 0; pos--)
				{
					if (item[pos-1] != ' ' && item[pos-1] != '\t') break;
				}
				item[pos] = '\0';
				
				if (strcmp(item, "Phi") == 0 || strcmp(item, "phi") == 0)
					pvalue |= MASK_PHI;
				else if (strcmp(item, "Chi") == 0 || strcmp(item, "chi") == 0)
					pvalue |= MASK_CHI;
				else if (strcmp(item, "Pot") == 0 || strcmp(item, "pot") == 0 || strcmp(item, "Psi_N") == 0 || strcmp(item, "psi_N") == 0 || strcmp(item, "PsiN") == 0 || strcmp(item, "psiN") == 0)
					pvalue |= MASK_POT;
				else if (strcmp(item, "B") == 0 || strcmp(item, "Bi") == 0)
					pvalue |= MASK_B;
				else if (strcmp(item, "P") == 0 || strcmp(item, "p") == 0 || strcmp(item, "v") == 0)
					pvalue |= MASK_P;
				else if (strcmp(item, "T00") == 0 || strcmp(item, "rho") == 0)
					pvalue |= MASK_T00;
				else if (strcmp(item, "Tij") == 0)
					pvalue |= MASK_TIJ;
				else if (strcmp(item, "rho_N") == 0 || strcmp(item, "rhoN") == 0)
					pvalue |= MASK_RBARE;
				else if (strcmp(item, "hij") == 0 || strcmp(item, "GW") == 0)
					pvalue |= MASK_HIJ;
				else if (strcmp(item, "Gadget") == 0 || strcmp(item, "Gadget2") == 0 || strcmp(item, "gadget") == 0 || strcmp(item, "gadget2") == 0)
					pvalue |= MASK_GADGET;
				else if (strcmp(item, "Particles") == 0 || strcmp(item, "particles") == 0 || strcmp(item, "pcls") == 0 || strcmp(item, "part") == 0)
					pvalue |= MASK_PCLS;
				else if (strcmp(item, "cross") == 0 || strcmp(item, "X-spectra") == 0 || strcmp(item, "x-spectra") == 0)
					pvalue |= MASK_XSPEC;
				else if (strcmp(item, "delta") == 0 || strcmp(item, "Ds") == 0 || strcmp(item, "D_s") == 0)
					pvalue |= MASK_DELTA;
				else if (strcmp(item, "delta_N") == 0 || strcmp(item, "deltaN") == 0)
					pvalue |= MASK_DBARE;
					
				start = comma+1;
				while (*start == ' ' || *start == '\t') start++;
			}  
			
			if (strcmp(start, "Phi") == 0 || strcmp(start, "phi") == 0)
				pvalue |= MASK_PHI;
			else if (strcmp(start, "Chi") == 0 || strcmp(start, "chi") == 0)
				pvalue |= MASK_CHI;
			else if (strcmp(start, "Pot") == 0 || strcmp(start, "pot") == 0 || strcmp(start, "Psi_N") == 0 || strcmp(start, "psi_N") == 0 || strcmp(start, "PsiN") == 0 || strcmp(start, "psiN") == 0)
				pvalue |= MASK_POT;
			else if (strcmp(start, "B") == 0 || strcmp(start, "Bi") == 0)
				pvalue |= MASK_B;
			else if (strcmp(start, "P") == 0 || strcmp(start, "p") == 0 || strcmp(start, "v") == 0)
				pvalue |= MASK_P;
			else if (strcmp(start, "T00") == 0 || strcmp(start, "rho") == 0)
				pvalue |= MASK_T00;
			else if (strcmp(start, "Tij") == 0)
				pvalue |= MASK_TIJ;
			else if (strcmp(start, "rho_N") == 0 || strcmp(start, "rhoN") == 0)
				pvalue |= MASK_RBARE;
			else if (strcmp(start, "hij") == 0 || strcmp(start, "GW") == 0)
				pvalue |= MASK_HIJ;
			else if (strcmp(start, "Gadget") == 0 || strcmp(start, "Gadget2") == 0 || strcmp(start, "gadget") == 0 || strcmp(start, "gadget2") == 0)
				pvalue |= MASK_GADGET;
			else if (strcmp(start, "Particles") == 0 || strcmp(start, "particles") == 0 || strcmp(start, "pcls") == 0 || strcmp(start, "part") == 0)
				pvalue |= MASK_PCLS;
			else if (strcmp(start, "cross") == 0 || strcmp(start, "X-spectra") == 0 || strcmp(start, "x-spectra") == 0)
				pvalue |= MASK_XSPEC;
			else if (strcmp(start, "delta") == 0 || strcmp(start, "Ds") == 0 || strcmp(start, "D_s") == 0)
				pvalue |= MASK_DELTA;
			else if (strcmp(start, "delta_N") == 0 || strcmp(start, "deltaN") == 0)
				pvalue |= MASK_DBARE;
			
			params[i].used = true;
			return true;
		}
	}
	
	return false;
}


//////////////////////////
// parseMetadata
//////////////////////////
// Description:
//   parses all metadata from the parameter array
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   sim        reference to metadata stucture (holds simulation parameters)
//   cosmo      reference to cosmology structure (holds cosmological parameters)
//   ic         reference to icsettings structure (holds settings for IC generation)
//
// Returns:
//   number of parameters parsed
// 
//////////////////////////

int parseMetadata(parameter * & params, const int numparam, metadata & sim, cosmology & cosmo, icsettings & ic)
{
	char par_string[PARAM_MAX_LENGTH];
	char * ptr;
	char * pptr[MAX_PCL_SPECIES];
	int usedparams = 0;
	int i;

	// parse settings for IC generator
	
	ic.pkfile[0] = '\0';
	ic.tkfile[0] = '\0';
	ic.metricfile[0][0] = '\0';
	ic.metricfile[1][0] = '\0';
	ic.metricfile[2][0] = '\0';
	ic.seed = 0;
	ic.flags = 0;
	ic.z_ic = -2.;
	ic.z_relax = -2.;
	ic.Cf = -1.0;
	ic.A_s = P_SPECTRAL_AMP;
	ic.n_s = P_SPECTRAL_INDEX;
	ic.k_pivot = P_PIVOT_SCALE;
	ic.restart_cycle = -1;
	ic.restart_tau = 0.;
	ic.restart_dtau = 0.;
	ic.restart_version = -1.;
	
	parseParameter(params, numparam, "seed", ic.seed);
	
	if (parseParameter(params, numparam, "IC generator", par_string))
	{
		if (par_string[0] == 'B' || par_string[0] == 'b')
			ic.generator = ICGEN_BASIC;
		else if (par_string[0] == 'R' || par_string[0] == 'r')
			ic.generator = ICGEN_READ_FROM_DISK;
#ifdef ICGEN_PREVOLUTION
		else if (par_string[0] == 'P' || par_string[0] == 'p')
			ic.generator = ICGEN_PREVOLUTION;
#endif
#ifdef ICGEN_SONG
		else if (par_string[0] == 'S' || par_string[0] == 's')
			ic.generator = ICGEN_SONG;
#endif
#ifdef ICGEN_FALCONIC
		else if (par_string[0] == 'F' || par_string[0] == 'f')
			ic.generator = ICGEN_FALCONIC;
#endif
		else
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": IC generator not recognized!" << endl;
			parallel.abortForce();
		}
	}
	else
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": IC generator not specified, selecting default (basic)" << endl;
		ic.generator = ICGEN_BASIC;
	}
	
	for (i = 0; i < MAX_PCL_SPECIES; i++)
		pptr[i] = ic.pclfile[i];
	if (ic.generator == ICGEN_READ_FROM_DISK)
	{
		if (!parseParameter(params, numparam, "particle file", pptr, i))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no particle file specified!" << endl;
			parallel.abortForce();
		}
	}
	else if (!parseParameter(params, numparam, "template file", pptr, i))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no template file specified!" << endl;
		parallel.abortForce();
	}
	
	for (; i < MAX_PCL_SPECIES; i++)
		strcpy(ic.pclfile[i], ic.pclfile[i-1]);

#ifdef ICGEN_FALCONIC
	if ((!parseParameter(params, numparam, "mPk file", ic.pkfile) && !parseParameter(params, numparam, "Tk file", ic.tkfile) && ic.generator != ICGEN_READ_FROM_DISK && ic.generator != ICGEN_FALCONIC)
#else
	if ((!parseParameter(params, numparam, "mPk file", ic.pkfile) && !parseParameter(params, numparam, "Tk file", ic.tkfile) && ic.generator != ICGEN_READ_FROM_DISK)
#endif
#ifdef ICGEN_PREVOLUTION
	    || ic.generator == ICGEN_PREVOLUTION)
#else
		)
#endif
	{
#ifdef HAVE_CLASS
		COUT << " initial transfer functions will be computed by calling CLASS" << endl;
#else
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no power spectrum file nor transfer function file specified!" << endl;
		parallel.abortForce();
#endif
	}
	
	if (parseParameter(params, numparam, "correct displacement", par_string))
	{
		if (par_string[0] == 'Y' || par_string[0] == 'y')
			ic.flags |= ICFLAG_CORRECT_DISPLACEMENT;
		else if (par_string[0] != 'N' && par_string[0] != 'n')
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": setting chosen for deconvolve displacement option not recognized, using default (no)" << endl;
	}
	
	if (parseParameter(params, numparam, "k-domain", par_string))
	{
		if (par_string[0] == 'S' || par_string[0] == 's')
			ic.flags |= ICFLAG_KSPHERE;
		else if (par_string[0] != 'C' && par_string[0] != 'c')
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": setting chosen for k-domain option not recognized, using default (cube)" << endl;
	}
	
	for (i = 0; i < MAX_PCL_SPECIES; i++)
		ic.numtile[i] = 0;
	
	if(!parseParameter(params, numparam, "tiling factor", ic.numtile, i) && ic.generator != ICGEN_READ_FROM_DISK)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tiling factor not specified, using default value for all species (1)" << endl;
		ic.numtile[0] = 1;
		i = 1;
	}
	
	for (; i < MAX_PCL_SPECIES; i++)
		ic.numtile[i] = ic.numtile[i-1];
	
	for (i = 0; i < MAX_PCL_SPECIES; i++)
	{
		if (ic.numtile[i] <= 0 && ic.generator != ICGEN_READ_FROM_DISK)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tiling number for particle template not set properly; using default value (1)" << endl;
			ic.numtile[i] = 1;
		}
	}
	
	if (ic.pkfile[0] != '\0')
	{
		sim.baryon_flag = 0;
		if (parseParameter(params, numparam, "baryon treatment", par_string))
		{
			if (par_string[0] != 'i' && par_string[0] != 'I')
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using mPk file forces baryon treatment = ignore" << endl;
			}
		}
	}
	else if (parseParameter(params, numparam, "baryon treatment", par_string))
	{
		if (par_string[0] == 'i' || par_string[0] == 'I')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "ignore" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 0;
		}
		else if (par_string[0] == 's' || par_string[0] == 'S')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "sample" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 1;
		}
		else if (par_string[0] == 'b' || par_string[0] == 'B')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "blend" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 2;
		}
		else if (par_string[0] == 'h' || par_string[0] == 'H')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "hybrid" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 3;
		}
		else
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": baryon treatment not supported!" << endl;
			parallel.abortForce();
		}
	}
	else if (ic.generator == ICGEN_READ_FROM_DISK)
	{
		sim.baryon_flag = 0;
	}
	else
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": baryon treatment not specified, using default (blend)" << endl;
		sim.baryon_flag = 2;
	}
	
	if (parseParameter(params, numparam, "radiation treatment", par_string))
	{
		if (par_string[0] == 'i' || par_string[0] == 'I')
		{
			sim.radiation_flag = 0;
			COUT << " radiation treatment set to: " << COLORTEXT_CYAN << "ignore" << COLORTEXT_RESET << endl;
		}
#ifdef HAVE_CLASS
		else if (par_string[0] == 'c' || par_string[0] == 'C')
		{
			sim.radiation_flag = 1;
			COUT << " radiation treatment set to: " << COLORTEXT_CYAN << "CLASS" << COLORTEXT_RESET << endl;
			if ((ic.pkfile[0] != '\0' || ic.tkfile[0] != '\0')
#ifdef ICGEN_PREVOLUTION
				&& ic.generator != ICGEN_PREVOLUTION
#endif
				)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using radiation treatment = CLASS and providing initial power spectra / transfer functions independently" << endl;
				COUT << "              is dangerous! In order to ensure consistency, it is recommended to call CLASS directly." << endl;
			}
		}
		else
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": radiation treatment not supported!" << endl;
			parallel.abortForce();
		}
#else
		else
		{
			sim.radiation_flag = 0;
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": CLASS is not available, setting radiation treatment = ignore" << endl;
		}
#endif
	}
	else
			sim.radiation_flag = 0;
	
	parseParameter(params, numparam, "relaxation redshift", ic.z_relax);

	if (ic.generator == ICGEN_READ_FROM_DISK)
	{
		parseParameter(params, numparam, "restart redshift", ic.z_ic);
		if(!parseParameter(params, numparam, "cycle", ic.restart_cycle))
		{
			if (sim.radiation_flag)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using radiation treatment = CLASS and IC generator = read from disk does not guarantee" << endl;
				COUT << "              that the realization of the radiation perturbations is consistent!" << endl;
			}
		}
		parseParameter(params, numparam, "tau", ic.restart_tau);
		parseParameter(params, numparam, "dtau", ic.restart_dtau);
		for (i = 0; i < 3; i++)
			pptr[i] = ic.metricfile[i];
		parseParameter(params, numparam, "metric file", pptr, i);
		if (parseParameter(params, numparam, "gevolution version", ic.restart_version))
		{
			if (ic.restart_version - GEVOLUTION_VERSION > 0.0001)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": version number of settings file (" << ic.restart_version << ") is higher than version of executable (" << GEVOLUTION_VERSION << ")!" << endl;
			}
		}
	}
#ifdef ICGEN_PREVOLUTION	
	else if (ic.generator == ICGEN_PREVOLUTION)
	{	
		if (!parseParameter(params, numparam, "prevolution redshift", ic.z_ic))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no starting redshift specified for IC generator = prevolution" << endl;
			parallel.abortForce();
		}
		
		parseParameter(params, numparam, "prevolution Courant factor", ic.Cf);
	}
#endif

	if (!parseParameter(params, numparam, "A_s", ic.A_s) && (
#ifdef ICGEN_FALCONIC
		ic.generator == ICGEN_FALCONIC ||
#endif
#ifdef ICGEN_PREVOLUTION
		ic.generator == ICGEN_PREVOLUTION ||
#endif		
		sim.radiation_flag > 0 || (ic.pkfile[0] == '\0' && ic.generator != ICGEN_READ_FROM_DISK)))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": power spectrum normalization not specified, using default value (2.215e-9)" << endl;
	}

	if (!parseParameter(params, numparam, "n_s", ic.n_s) && (	
#ifdef ICGEN_FALCONIC
		ic.generator == ICGEN_FALCONIC ||
#endif
#ifdef ICGEN_PREVOLUTION
		ic.generator == ICGEN_PREVOLUTION ||
#endif
		sim.radiation_flag > 0 || (ic.pkfile[0] == '\0' && ic.generator != ICGEN_READ_FROM_DISK)))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": scalar spectral index not specified, using default value (0.9619)" << endl;
	}
	
	if (!parseParameter(params, numparam, "k_pivot", ic.k_pivot) && (
#ifdef ICGEN_FALCONIC
		ic.generator == ICGEN_FALCONIC ||
#endif
#ifdef ICGEN_PREVOLUTION
		ic.generator == ICGEN_PREVOLUTION ||
#endif
		sim.radiation_flag > 0 || (ic.pkfile[0] == '\0' && ic.generator != ICGEN_READ_FROM_DISK)))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": pivot scale not specified, using default value (0.05 / Mpc)" << endl;
	}
	
	// parse metadata
	
	sim.numpts = 0;
	sim.downgrade_factor = 1;
	for (i = 0; i < MAX_PCL_SPECIES; i++) sim.numpcl[i] = 0;
	sim.vector_flag = VECTOR_PARABOLIC;
	sim.gr_flag = 0;
	sim.out_pk = 0;
	sim.out_snapshot = 0;
	sim.num_pk = MAX_OUTPUTS;
	sim.numbins = 0;
	sim.num_snapshot = MAX_OUTPUTS;
	sim.num_restart = MAX_OUTPUTS;
	for (i = 0; i < MAX_PCL_SPECIES; i++) sim.tracer_factor[i] = 1;
	sim.Cf = 1.;
	sim.steplimit = 1.;
	sim.boxsize = -1.;
	sim.wallclocklimit = -1.;
	sim.z_in = 0.;

	if (parseParameter(params, numparam, "vector method", par_string))
	{
		if (par_string[0] == 'p' || par_string[0] == 'P')
		{
			COUT << " vector method set to: " << COLORTEXT_CYAN << "parabolic" << COLORTEXT_RESET << endl;
			sim.vector_flag = VECTOR_PARABOLIC;
		}
		else if (par_string[0] == 'e' || par_string[0] == 'E')
		{
			COUT << " vector method set to: " << COLORTEXT_CYAN << "elliptic" << COLORTEXT_RESET << endl;
			sim.vector_flag = VECTOR_ELLIPTIC;
		}
		else
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": vector method not supported!" << endl;
			parallel.abortForce();
		}
	}
	
	if (!parseParameter(params, numparam, "generic file base", sim.basename_generic))
		sim.basename_generic[0] = '\0';
	
	if (!parseParameter(params, numparam, "snapshot file base", sim.basename_snapshot))
		strcpy(sim.basename_snapshot, "snapshot");
		
	if (!parseParameter(params, numparam, "Pk file base", sim.basename_pk))
		strcpy(sim.basename_pk, "pk");
		
	if (!parseParameter(params, numparam, "output path", sim.output_path))
		sim.output_path[0] = '\0';
		
	if (!parseParameter(params, numparam, "hibernation path", sim.restart_path))
		strcpy(sim.restart_path, sim.output_path);
		
	if (!parseParameter(params, numparam, "hibernation file base", sim.basename_restart))
		strcpy(sim.basename_restart, "restart");
		
	parseParameter(params, numparam, "boxsize", sim.boxsize);
	if (sim.boxsize <= 0. || !isfinite(sim.boxsize))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": simulation box size not set properly!" << endl;
		parallel.abortForce();
	}
	
	parseParameter(params, numparam, "Ngrid", sim.numpts);
	if (sim.numpts < 2 || !isfinite(sim.numpts))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of grid points not set properly!" << endl;
		parallel.abortForce();
	}

	if (parseParameter(params, numparam, "downgrade factor", sim.downgrade_factor))
	{
		if (sim.downgrade_factor < 1 || sim.downgrade_factor >= sim.numpts || !isfinite(sim.downgrade_factor))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": downgrade factor makes no sense!" << endl;
			parallel.abortForce();
		}
		if (sim.downgrade_factor > 1 && (sim.numpts % parallel.grid_size()[0] || sim.numpts % parallel.grid_size()[1] || (sim.numpts / parallel.grid_size()[0]) % sim.downgrade_factor || (sim.numpts / parallel.grid_size()[1]) % sim.downgrade_factor))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": downgrade factor does not appear to be compatible with process layout at given Ngrid!" << endl;
			parallel.abortForce();
		}
	}
	
	if (!parseParameter(params, numparam, "move limit", sim.movelimit))
		sim.movelimit = (double) sim.numpts;
	
	if (!parseParameter(params, numparam, "initial redshift", sim.z_in))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": initial redshift not specified!" << endl;
		parallel.abortForce();
	}
	
	if (ic.z_relax < -1.) ic.z_relax = sim.z_in;
#ifdef ICGEN_PREVOLUTION	
	else if (ic.generator == ICGEN_PREVOLUTION && ic.z_relax < sim.z_in)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": relaxation redshift cannot be below initial redshift for IC generator = prevolution; reset to initial redshift!" << endl;
		ic.z_relax = sim.z_in;
	}
#endif

	if (ic.z_ic < sim.z_in && ic.generator != ICGEN_READ_FROM_DISK) ic.z_ic = sim.z_in;
	
	parseParameter(params, numparam, "snapshot redshifts", sim.z_snapshot, sim.num_snapshot);
	if (sim.num_snapshot > 0)
		qsort((void *) sim.z_snapshot, (size_t) sim.num_snapshot, sizeof(double), sort_descending);
	
	parseParameter(params, numparam, "Pk redshifts", sim.z_pk, sim.num_pk);
	if (sim.num_pk > 0)
		qsort((void *) sim.z_pk, (size_t) sim.num_pk, sizeof(double), sort_descending);
		
	parseParameter(params, numparam, "hibernation redshifts", sim.z_restart, sim.num_restart);
	if (sim.num_restart > 0)
		qsort((void *) sim.z_restart, (size_t) sim.num_restart, sizeof(double), sort_descending);
		
	parseParameter(params, numparam, "hibernation wallclock limit", sim.wallclocklimit);
		
	parseFieldSpecifiers(params, numparam, "snapshot outputs", sim.out_snapshot);
	parseFieldSpecifiers(params, numparam, "Pk outputs", sim.out_pk);
	
	i = MAX_PCL_SPECIES;
	parseParameter(params, numparam, "tracer factor", sim.tracer_factor, i);
	for (; i > 0; i--)
	{
		if (sim.tracer_factor[i-1] < 1)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tracer factor not set properly; using default value (1)" << endl;
			sim.tracer_factor[i-1] = 1;
		}
	}
	
	if ((sim.num_snapshot <= 0 || sim.out_snapshot == 0) && (sim.num_pk <= 0 || sim.out_pk == 0))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": no output specified!" << endl;
	}
	
	if (!parseParameter(params, numparam, "Pk bins", sim.numbins))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of Pk bins not set properly; using default value (64)" << endl;
		sim.numbins = 64;
	}
	
	parseParameter(params, numparam, "Courant factor", sim.Cf);

	if (ic.Cf < 0.) ic.Cf = sim.Cf;
	
	parseParameter(params, numparam, "time step limit", sim.steplimit);
	
	if (parseParameter(params, numparam, "gravity theory", par_string))
	{
		if (par_string[0] == 'N' || par_string[0] == 'n')
		{
			COUT << " gravity theory set to: " << COLORTEXT_CYAN << "Newtonian" << COLORTEXT_RESET << endl;
			sim.gr_flag = 0;
			if (ic.pkfile[0] == '\0' && ic.tkfile[0] != '\0'
#ifdef ICGEN_PREVOLUTION
				&& ic.generator != ICGEN_PREVOLUTION
#endif
				)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": gauge transformation to N-body gauge can only be performed for the positions; the transformation for" << endl;
				COUT << "              the velocities requires time derivatives of transfer functions. Call CLASS directly to avoid this issue." << endl;
			}
		}
		else if (par_string[0] == 'G' || par_string[0] == 'g')
		{
			COUT << " gravity theory set to: " << COLORTEXT_CYAN << "General Relativity" << COLORTEXT_RESET << endl;
			sim.gr_flag = 1;
		}
		else
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": gravity theory unknown, using default (General Relativity)" << endl;
			sim.gr_flag = 1;
		}
	}
	else
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": gravity theory not selected, using default (General Relativity)" << endl;
		sim.gr_flag = 1;
	}
	
	
	// parse cosmological parameters
	
	if (!parseParameter(params, numparam, "h", cosmo.h))
	{
		cosmo.h = P_HUBBLE;
	}
	
	cosmo.num_ncdm = MAX_PCL_SPECIES-2;
	if (!parseParameter(params, numparam, "m_ncdm", cosmo.m_ncdm, cosmo.num_ncdm))
	{
		for (i = 0; i < MAX_PCL_SPECIES-2; i++) cosmo.m_ncdm[i] = 0.;
		cosmo.num_ncdm = 0;
	}
	
	if (parseParameter(params, numparam, "N_ncdm", i))
	{
		if (i < 0 || !isfinite(i))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of ncdm species not set properly!" << endl;
			parallel.abortForce();
		}
		if (i > cosmo.num_ncdm)
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": N_ncdm = " << i << " is larger than the number of mass parameters specified (" << cosmo.num_ncdm << ")!" << endl;
			parallel.abortForce();
		}
		cosmo.num_ncdm = i;
	}
	else if (cosmo.num_ncdm > 0)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": N_ncdm not specified, inferring from number of mass parameters in m_ncdm (" << cosmo.num_ncdm << ")!" << endl;
	}
	
	for (i = 0; i < MAX_PCL_SPECIES-2; i++)
	{
		cosmo.T_ncdm[i] = P_T_NCDM;
		cosmo.deg_ncdm[i] = 1.0;	
	}
	parseParameter(params, numparam, "T_ncdm", cosmo.T_ncdm, i);
	i = MAX_PCL_SPECIES-2;
	parseParameter(params, numparam, "deg_ncdm", cosmo.deg_ncdm, i);
	
	for (i = 0; i < cosmo.num_ncdm; i++)
	{
		cosmo.Omega_ncdm[i] = cosmo.m_ncdm[i] * cosmo.deg_ncdm[i] / P_NCDM_MASS_OMEGA / cosmo.h / cosmo.h;
	}
	
	if (parseParameter(params, numparam, "T_cmb", cosmo.Omega_g))
	{
		cosmo.Omega_g = cosmo.Omega_g * cosmo.Omega_g / cosmo.h;
		cosmo.Omega_g = cosmo.Omega_g * cosmo.Omega_g * C_PLANCK_LAW; // Planck's law
	}
	else if (parseParameter(params, numparam, "omega_g", cosmo.Omega_g))
	{
		cosmo.Omega_g /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_g", cosmo.Omega_g))
	{
		cosmo.Omega_g = 0.;
	}
	
	if (parseParameter(params, numparam, "N_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur *= (7./8.) * pow(4./11., 4./3.) * cosmo.Omega_g;
	}
	else if (parseParameter(params, numparam, "N_eff", cosmo.Omega_ur))
	{
		cosmo.Omega_ur *= (7./8.) * pow(4./11., 4./3.) * cosmo.Omega_g;
	}
	else if (parseParameter(params, numparam, "omega_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur = P_N_UR * (7./8.) * pow(4./11., 4./3.) * cosmo.Omega_g;
	}
	
	cosmo.Omega_rad = cosmo.Omega_g + cosmo.Omega_ur;
	
	if (parseParameter(params, numparam, "omega_b", cosmo.Omega_b))
	{
		cosmo.Omega_b /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_b", cosmo.Omega_b))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": Omega_b not found in settings file, setting to default (0)." << endl;
		cosmo.Omega_b = 0.;
	}
	
	if (parseParameter(params, numparam, "omega_cdm", cosmo.Omega_cdm))
	{
		cosmo.Omega_cdm /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_cdm", cosmo.Omega_cdm))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": Omega_cdm not found in settings file, setting to default (1)." << endl;
		cosmo.Omega_cdm = 1.;
	}
	
	cosmo.Omega_m = cosmo.Omega_cdm + cosmo.Omega_b;
	for (i = 0; i < cosmo.num_ncdm; i++) cosmo.Omega_m += cosmo.Omega_ncdm[i];
	
	if (cosmo.Omega_m <= 0. || cosmo.Omega_m > 1.)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": total matter density out of range!" << endl;
		parallel.abortForce();
	}
	else if (cosmo.Omega_rad < 0. || cosmo.Omega_rad > 1. - cosmo.Omega_m)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": total radiation energy density out of range!" << endl;
		parallel.abortForce();
	}
	else
	{
		COUT << " cosmological parameters are: Omega_m0 = " << cosmo.Omega_m << ", Omega_rad0 = " << cosmo.Omega_rad << ", h = " << cosmo.h << endl;
		cosmo.Omega_Lambda = 1. - cosmo.Omega_m - cosmo.Omega_rad;
	}

	if(!parseParameter(params, numparam, "switch delta_rad", sim.z_switch_deltarad))
		sim.z_switch_deltarad = 0.;

	i = MAX_PCL_SPECIES-2;
	if (!parseParameter(params, numparam, "switch delta_ncdm", sim.z_switch_deltancdm, i))
	{
		sim.z_switch_deltancdm[0] = sim.z_in;
		i = 1;
	}
	for (; i < MAX_PCL_SPECIES-2; i++)
		sim.z_switch_deltancdm[i] = sim.z_switch_deltancdm[i-1];

	if(!parseParameter(params, numparam, "switch linear chi", sim.z_switch_linearchi))
	{
		if (sim.gr_flag > 0)
			sim.z_switch_linearchi = 0.;
		else
			sim.z_switch_linearchi = 0.011;
			
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.z_switch_linearchi < sim.z_switch_deltancdm[i])
				sim.z_switch_linearchi = sim.z_switch_deltancdm[i];
		}
	}
	else if (sim.gr_flag == 0 && sim.z_switch_linearchi <= 0.01)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": with garavity theory = Newton the switch linear chi redshift must be larger than 0.01." << endl;
		COUT << "              setting switch linear chi = 0.011" << endl;
		sim.z_switch_linearchi = 0.011;
	}

	i = MAX_PCL_SPECIES-2;
	if (!parseParameter(params, numparam, "switch B ncdm", sim.z_switch_Bncdm, i))
	{
		for (i = 0; i < MAX_PCL_SPECIES-2; i++)
			sim.z_switch_Bncdm[i] = sim.z_switch_deltancdm[i];
	}
	for (; i < MAX_PCL_SPECIES-2; i++)
		sim.z_switch_Bncdm[i] = sim.z_switch_Bncdm[i-1];
	
	for (i = 0; i < numparam; i++)
	{
		if (params[i].used) usedparams++;
	}
	
	return usedparams;
}

#endif
