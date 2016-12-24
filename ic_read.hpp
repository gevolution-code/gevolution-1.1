//////////////////////////
// ic_read.hpp
//////////////////////////
// 
// read initial conditions from disk
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: December 2016
//
//////////////////////////

#ifndef IC_READ_HEADER
#define IC_READ_HEADER

//////////////////////////
// readIC
//////////////////////////
// Description:
//   reads initial conditions from disk
//
// Arguments: 
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              reference to scale factor
//   tau            reference to conformal coordinate time
//   dtau           time step
//   dtau_old       previous time step (will be passed back)
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
//   maxvel         array that will contain the maximum q/m/a (max. velocity)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   cycle          reference to cycle counter (for restart after hibernation)
//   snapcount      reference to snapshot counter (for restart after hibernation)
//   pkcount        reference to spectra counter (for restart after hibernation)
//   restartcount   reference to restart counter (for restart after hibernation)
//
// Returns:
// 
//////////////////////////

void readIC(metadata & sim, icsettings & ic, cosmology & cosmo, const double fourpiG, double & a, double & tau, double & dtau, double & dtau_old, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, double * maxvel, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, int & cycle, int & snapcount, int & pkcount, int & restartcount)
{
	part_simple_info pcls_cdm_info;
	part_simple_dataType pcls_cdm_dataType;
	part_simple_info pcls_b_info;
	part_simple_dataType pcls_b_dataType;
	part_simple_info pcls_ncdm_info[MAX_PCL_SPECIES];
	part_simple_dataType pcls_ncdm_dataType;
	Real boxSize[3] = {1.,1.,1.};
	string filename;
	string buf;
	int i, p;
	char * ext;
	char line[PARAM_MAX_LINESIZE];
	FILE * bgfile;
	struct fileDsc fd;
	gadget2_header hdr;
	long * numpcl;
	Real * dummy1;
	Real * dummy2;
	Site x(Bi->lattice());
	rKSite kFT(scalarFT->lattice());
	
	filename.reserve(PARAM_MAX_LENGTH);
	hdr.npart[1] = 0;
	
	projection_init(phi);
	
	if (ic.z_ic > -1.)
		a = 1. / (1. + ic.z_ic);
	
	strcpy(pcls_cdm_info.type_name, "part_simple");
	pcls_cdm_info.mass = 0.;
	pcls_cdm_info.relativistic = false;
	
	pcls_cdm->initialize(pcls_cdm_info, pcls_cdm_dataType, &(phi->lattice()), boxSize);
	
	if ((ext = strstr(ic.pclfile[0], ".h5")) != NULL)
	{
		filename.assign(ic.pclfile[0], ext-ic.pclfile[0]);
		get_fileDsc_global(filename + ".h5", fd);
		numpcl = (long *) malloc(fd.numProcPerFile * sizeof(long));
		dummy1 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
		dummy2 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
		get_fileDsc_local(filename + ".h5", numpcl, dummy1, dummy2, fd.numProcPerFile);
		for (i = 0; i < fd.numProcPerFile; i++)
			sim.numpcl[0] += numpcl[i];
		pcls_cdm->loadHDF5(filename, 1);
		free(numpcl);
		free(dummy1);
		free(dummy2);
	}
	else
	{
		i = 0;
		do
		{
			filename.assign(ic.pclfile[0]);
			pcls_cdm->loadGadget2(filename, hdr);
			if (hdr.npart[1] == 0) break;
			if (hdr.time / a > 1.001 || hdr.time / a < 0.999)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": redshift indicated in Gadget2 header does not match initial redshift of simulation!" << endl;
			}
			sim.numpcl[0] += hdr.npart[1];
			i++;
			if (hdr.num_files > 1)
			{
				ext = ic.pclfile[0];
				while (strchr(ext, (int) '.') != NULL)
					ext = strchr(ext, (int) '.');
				sprintf(ext+1, "%d", i);
			}
		}
		while (i < hdr.num_files);
		
		if (sim.baryon_flag == 1)
			pcls_cdm->parts_info()->mass = cosmo.Omega_cdm / (Real) sim.numpcl[0];
		else
			pcls_cdm->parts_info()->mass = (cosmo.Omega_cdm + cosmo.Omega_b) / (Real) sim.numpcl[0];
	}
	
	COUT << " " << sim.numpcl[0] << " cdm particles read successfully." << endl;
	maxvel[0] = pcls_cdm->updateVel(update_q, 0., &phi, 1, &a);
	
	if (sim.baryon_flag == 1)
	{
		strcpy(pcls_b_info.type_name, "part_simple");
		pcls_b_info.mass = 0.;
		pcls_b_info.relativistic = false;
	
		pcls_b->initialize(pcls_b_info, pcls_b_dataType, &(phi->lattice()), boxSize);
		
		if ((ext = strstr(ic.pclfile[1], ".h5")) != NULL)
		{
			filename.assign(ic.pclfile[1], ext-ic.pclfile[1]);
			get_fileDsc_global(filename + ".h5", fd);
			numpcl = (long *) malloc(fd.numProcPerFile * sizeof(long));
			dummy1 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			dummy2 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			get_fileDsc_local(filename + ".h5", numpcl, dummy1, dummy2, fd.numProcPerFile);
			for (i = 0; i < fd.numProcPerFile; i++)
				sim.numpcl[1] += numpcl[i];
			pcls_b->loadHDF5(filename, 1);
			free(numpcl);
			free(dummy1);
			free(dummy2);
		}
		else
		{
			i = 0;
			do
			{
				filename.assign(ic.pclfile[1]);
				pcls_b->loadGadget2(filename, hdr);
				if (hdr.npart[1] == 0) break;
				if (hdr.time / a > 1.001 || hdr.time / a < 0.999)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": redshift indicated in Gadget2 header does not match initial redshift of simulation!" << endl;
				}
				sim.numpcl[1] += hdr.npart[1];
				i++;
				if (hdr.num_files > 1)
				{
					ext = ic.pclfile[1];
					while (strchr(ext, (int) '.') != NULL)
						ext = strchr(ext, (int) '.');
					sprintf(ext+1, "%d", i);
				}
			}
			while (i < hdr.num_files);
			
			pcls_b->parts_info()->mass = cosmo.Omega_b / (Real) sim.numpcl[1];
		}
		
		COUT << " " << sim.numpcl[1] << " baryon particles read successfully." << endl;
		maxvel[1] = pcls_b->updateVel(update_q, 0., &phi, 1, &a);
	}
	else
		sim.baryon_flag = 0;
		
	for (p = 0; p < cosmo.num_ncdm; p++)
	{
		strcpy(pcls_ncdm_info[p].type_name, "part_simple");
		pcls_ncdm_info[p].mass = 0.;
		pcls_ncdm_info[p].relativistic = true;
	
		pcls_ncdm[p].initialize(pcls_ncdm_info[p], pcls_ncdm_dataType, &(phi->lattice()), boxSize);
		
		if ((ext = strstr(ic.pclfile[sim.baryon_flag+1+p], ".h5")) != NULL)
		{
			filename.assign(ic.pclfile[sim.baryon_flag+1+p], ext-ic.pclfile[sim.baryon_flag+1+p]);
			get_fileDsc_global(filename + ".h5", fd);
			numpcl = (long *) malloc(fd.numProcPerFile * sizeof(long));
			dummy1 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			dummy2 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			get_fileDsc_local(filename + ".h5", numpcl, dummy1, dummy2, fd.numProcPerFile);
			for (i = 0; i < fd.numProcPerFile; i++)
				sim.numpcl[1] += numpcl[i];
			pcls_ncdm[p].loadHDF5(filename, 1);
			free(numpcl);
			free(dummy1);
			free(dummy2);
		}
		else
		{
			i = 0;
			do
			{
				filename.assign(ic.pclfile[sim.baryon_flag+1+p]);
				pcls_ncdm[p].loadGadget2(filename, hdr);
				if (hdr.npart[1] == 0) break;
				if (hdr.time / a > 1.001 || hdr.time / a < 0.999)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": redshift indicated in Gadget2 header does not match initial redshift of simulation!" << endl;
				}
				sim.numpcl[sim.baryon_flag+1+p] += hdr.npart[1];
				i++;
				if (hdr.num_files > 1)
				{
					ext = ic.pclfile[sim.baryon_flag+1+p];
					while (strchr(ext, (int) '.') != NULL)
						ext = strchr(ext, (int) '.');
					sprintf(ext+1, "%d", i);
				}
			}
			while (i < hdr.num_files);
			
			pcls_ncdm[p].parts_info()->mass = cosmo.Omega_ncdm[p] / (Real) sim.numpcl[sim.baryon_flag+1+p];
		}
		
		COUT << " " << sim.numpcl[sim.baryon_flag+1+p] << " ncdm particles read successfully." << endl;
		maxvel[sim.baryon_flag+1+p] = pcls_ncdm[p].updateVel(update_q, 0., &phi, 1, &a);
	}
	
	if (sim.gr_flag > 0 && ic.metricfile[0][0] != '\0')
	{
		filename.assign(ic.metricfile[0]);
		phi->loadHDF5(filename);
	}
	else
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if (sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		scalarProjectionCIC_comm(source);
	
		plan_source->execute(FFT_FORWARD);
	
		kFT.first();
		if (kFT.coord(0) == 0 && kFT.coord(1) == 0 && kFT.coord(2) == 0)
			(*scalarFT)(kFT) = Cplx(0.,0.);
				
		solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a, 3. * sim.gr_flag * (Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) + fourpiG * cosmo.Omega_m / a));
		plan_phi->execute(FFT_BACKWARD);
	}

	phi->updateHalo();

	if (ic.restart_cycle >= 0)
	{	
#ifndef CHECK_B	
		if (sim.vector_flag == VECTOR_PARABOLIC)
#endif
		{
			filename.assign(ic.metricfile[2*sim.gr_flag]);
			Bi->loadHDF5(filename);

			for (x.first(); x.test(); x.next())
			{
				(*Bi)(x,0) *= a * a / (sim.numpts * sim.numpts);
				(*Bi)(x,1) *= a * a / (sim.numpts * sim.numpts);
				(*Bi)(x,2) *= a * a / (sim.numpts * sim.numpts);
			}
			plan_Bi->execute(FFT_FORWARD);
		}
		
		if (sim.gr_flag > 0)
		{
			filename.assign(ic.metricfile[1]);
			chi->loadHDF5(filename);
			chi->updateHalo();
		}
		
		if (parallel.isRoot())
		{
			sprintf(line, "%s%s_background.dat", sim.output_path, sim.basename_generic);
			bgfile = fopen(line, "r");
			if (bgfile == NULL)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to locate file for background output! A new file will be created" << endl;
				bgfile = fopen(line, "w");
				if (bgfile == NULL)
				{
					COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to create file for background output!" << endl;
					parallel.abortForce();
				}
				else
				{
					fprintf(bgfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0  phi(k=0)       T00(k=0)\n");
					fclose(bgfile);
				}
			}
			else
			{
				buf.reserve(PARAM_MAX_LINESIZE);
				buf.clear();
				
				if (fgets(line, PARAM_MAX_LINESIZE, bgfile) == 0)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read file for background output! A new file will be created" << endl;
				}
				else if (line[0] != '#')
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": file for background output has unexpected format! Contents will be overwritten!" << endl;
				}
				else
				{
					if (fgets(line, PARAM_MAX_LINESIZE, bgfile) == 0)
					{
						COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read file for background output! A new file will be created" << endl;
					}
					else if (line[0] != '#')
					{
						COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": file for background output has unexpected format! Contents will be overwritten!" << endl;
					}
				}
				
				while (fgets(line, PARAM_MAX_LINESIZE, bgfile) != 0)
				{
					if (sscanf(line, " %d", &i) != 1)
						break;
					
					if (i > ic.restart_cycle)
						break;
					else
						buf += line;
				}
	
				fclose(bgfile);
				
				sprintf(line, "%s%s_background.dat", sim.output_path, sim.basename_generic);
				bgfile = fopen(line, "w");
				
				if (bgfile == NULL)
				{
					COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to create file for background output!" << endl;
					parallel.abortForce();
				}
				else
				{
					fprintf(bgfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0  phi(k=0)       T00(k=0)\n");
					fwrite((const void *) buf.data(), sizeof(char), buf.length(), bgfile);
					fclose(bgfile);
					buf.clear();
				}
			}
		}
	}
	else
	{
		projection_init(Bi);
		projection_T0i_project(pcls_cdm, Bi, phi);
		if (sim.baryon_flag)
			projection_T0i_project(pcls_b, Bi, phi);
		projection_T0i_comm(Bi);
		plan_Bi->execute(FFT_FORWARD);
		projectFTvector(*BiFT, *BiFT, fourpiG / (double) sim.numpts / (double) sim.numpts);	
		plan_Bi->execute(FFT_BACKWARD);	
		Bi->updateHalo();
		
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if (sim.baryon_flag)
			projection_Tij_project(pcls_b, Sij, a, phi);
		projection_Tij_comm(Sij);
	
		prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / a / (double) sim.numpts / (double) sim.numpts);	
		plan_Sij->execute(FFT_FORWARD);	
		projectFTscalar(*SijFT, *scalarFT);
		plan_chi->execute(FFT_BACKWARD);		
		chi->updateHalo();
	}
		
	if (ic.restart_tau > 0.)
		tau = ic.restart_tau;
	else
		particleHorizon(a, fourpiG, cosmo);
		
	if (ic.restart_dtau > 0.)
		dtau_old = ic.restart_dtau;
		
	if (sim.Cf / (double) sim.numpts < sim.steplimit / Hconf(a, fourpiG, cosmo))
		dtau = sim.Cf / (double) sim.numpts;
	else
		dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);

	if (ic.restart_cycle >= 0)
		cycle = ic.restart_cycle + 1;
		
	while (snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
		snapcount++;
		
	while (pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
		pkcount++;
		
	while (restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
		restartcount++;
}

#endif

