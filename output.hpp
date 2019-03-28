//////////////////////////
// output.hpp
//////////////////////////
// 
// Output of snapshots and spectra
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: March 2019
//
//////////////////////////

#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using namespace std;


//////////////////////////
// writeSnapshots
//////////////////////////
// Description:
//   output of snapshots
// 
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   hdr            Gadget-2 header structure
//   a              scale factor
//   snapcount      snapshot index
//   h5filename     base name for HDF5 output file
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
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
//   Bi_check       pointer to allocated field (or NULL)
//   BiFT_check     pointer to allocated field (or NULL)
//   plan_Bi_check  pointer to FFT planner (or NULL)
//
// Returns:
// 
//////////////////////////

void writeSnapshots(metadata & sim, cosmology & cosmo, const double fourpiG, gadget2_header & hdr, const double a, const int snapcount, string h5filename, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, Field<Real> * Bi_check = NULL, Field<Cplx> * BiFT_check = NULL, PlanFFT<Cplx> * plan_Bi_check = NULL)
{
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[64];
	int i;
	Site x(phi->lattice());
	Real divB, curlB, divh, traceh, normh;

	sprintf(filename, "%03d", snapcount);
			
#ifdef EXTERNAL_IO
	while (ioserver.openOstream()== OSTREAM_FAIL);
	
	if (sim.out_snapshot & MASK_PCLS)
	{
		pcls_cdm->saveHDF5_server_open(h5filename + filename + "_cdm");
		if (sim.baryon_flag)
			pcls_b->saveHDF5_server_open(h5filename + filename + "_b");
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			pcls_ncdm[i].saveHDF5_server_open(h5filename + filename + buffer);
		}
	}
	
	if (sim.out_snapshot & MASK_T00)
		source->saveHDF5_server_open(h5filename + filename + "_T00");
				
	if (sim.out_snapshot & MASK_B)
		Bi->saveHDF5_server_open(h5filename + filename + "_B");
	
	if (sim.out_snapshot & MASK_PHI)
		phi->saveHDF5_server_open(h5filename + filename + "_phi");
				
	if (sim.out_snapshot & MASK_CHI)
		chi->saveHDF5_server_open(h5filename + filename + "_chi");
	
	if (sim.out_snapshot & MASK_HIJ)
		Sij->saveHDF5_server_open(h5filename + filename + "_hij");
				
#ifdef CHECK_B
	if (sim.out_snapshot & MASK_B)
		Bi_check->saveHDF5_server_open(h5filename + filename + "_B_check");
#endif
#endif		
			
	if (sim.out_snapshot & MASK_RBARE || sim.out_snapshot & MASK_POT)
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if (sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		for (i = 0; i < cosmo.num_ncdm; i++)
			scalarProjectionCIC_project(pcls_ncdm+i, source);
		scalarProjectionCIC_comm(source);
	}

	if (sim.out_snapshot & MASK_RBARE)
	{
		if (sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_rhoN.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_rhoN.h5");
	}
			
	if (sim.out_snapshot & MASK_POT)
	{
		plan_source->execute(FFT_FORWARD);				
		solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a);
		plan_source->execute(FFT_BACKWARD);
		if (sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_psiN.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_psiN.h5");
	}
				
	if (sim.out_snapshot & MASK_T00)
	{
		projection_init(source);
		if (sim.gr_flag > 0)
		{
			projection_T00_project(pcls_cdm, source, a, phi);
			if (sim.baryon_flag)
				projection_T00_project(pcls_b, source, a, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_T00_project(pcls_ncdm+i, source, a, phi);
		}
		else
		{
			scalarProjectionCIC_project(pcls_cdm, source);
			if (sim.baryon_flag)
				scalarProjectionCIC_project(pcls_b, source);
			for (i = 0; i < cosmo.num_ncdm; i++)
				scalarProjectionCIC_project(pcls_ncdm+i, source);
		}
		projection_T00_comm(source);
#ifdef EXTERNAL_IO
		source->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_T00.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_T00.h5");
#endif
	}
				
	if (sim.out_snapshot & MASK_B)
	{
		if (sim.gr_flag == 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
		}
		for (x.first(); x.test(); x.next())
		{
			(*Bi)(x,0) /= a * a * sim.numpts;
			(*Bi)(x,1) /= a * a * sim.numpts;
			(*Bi)(x,2) /= a * a * sim.numpts;
		}
		Bi->updateHalo();
				
		computeVectorDiagnostics(*Bi, divB, curlB);			
		COUT << " B diagnostics: max |divB| = " << divB << ", max |curlB| = " << curlB << endl;

#ifdef EXTERNAL_IO
		Bi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			Bi->saveHDF5_coarseGrain3D(h5filename + filename + "_B.h5", sim.downgrade_factor);
		else				
			Bi->saveHDF5(h5filename + filename + "_B.h5");
#endif
				
		if (sim.gr_flag > 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
			Bi->updateHalo();
		}
	}
			
	if (sim.out_snapshot & MASK_PHI)
#ifdef EXTERNAL_IO
		phi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			phi->saveHDF5_coarseGrain3D(h5filename + filename + "_phi.h5", sim.downgrade_factor);
		else
			phi->saveHDF5(h5filename + filename + "_phi.h5");
#endif
				
	if (sim.out_snapshot & MASK_CHI)
#ifdef EXTERNAL_IO
		chi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else	
		if (sim.downgrade_factor > 1)
			chi->saveHDF5_coarseGrain3D(h5filename + filename + "_chi.h5", sim.downgrade_factor);
		else
			chi->saveHDF5(h5filename + filename + "_chi.h5");
#endif
				
	if (sim.out_snapshot & MASK_HIJ)
	{
		projectFTtensor(*SijFT, *SijFT);
		plan_Sij->execute(FFT_BACKWARD);
		Sij->updateHalo();
				
		computeTensorDiagnostics(*Sij, divh, traceh, normh);
		COUT << " GW diagnostics: max |divh| = " << divh << ", max |traceh| = " << traceh << ", max |h| = " << normh << endl;

#ifdef EXTERNAL_IO
		Sij->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else	
		if (sim.downgrade_factor > 1)
			Sij->saveHDF5_coarseGrain3D(h5filename + filename + "_hij.h5", sim.downgrade_factor);
		else
			Sij->saveHDF5(h5filename + filename + "_hij.h5");
#endif
	}

	if (sim.out_snapshot & MASK_TIJ)
	{						
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if (sim.baryon_flag)
			projection_Tij_project(pcls_b, Sij, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_Tij_project(pcls_ncdm+i, Sij, a, phi);
		projection_Tij_comm(Sij);

		if (sim.downgrade_factor > 1)
			Sij->saveHDF5_coarseGrain3D(h5filename + filename + "_Tij.h5", sim.downgrade_factor);
		else
			Sij->saveHDF5(h5filename + filename + "_Tij.h5");
	}
			
	if (sim.out_snapshot & MASK_P)
	{
		projection_init(Bi);
		projection_T0i_project(pcls_cdm, Bi, phi);
		if (sim.baryon_flag)
			projection_T0i_project(pcls_b, Bi, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_T0i_project(pcls_ncdm+i, Bi, phi);
		projection_T0i_comm(Bi);
		if (sim.downgrade_factor > 1)
			Bi->saveHDF5_coarseGrain3D(h5filename + filename + "_p.h5", sim.downgrade_factor);
		else
			Bi->saveHDF5(h5filename + filename + "_p.h5");
		if (sim.gr_flag > 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
			Bi->updateHalo();
		}
	}
				
#ifdef CHECK_B
	if (sim.out_snapshot & MASK_B)
	{
		if (sim.vector_flag == VECTOR_PARABOLIC)
		{
			projection_init(Bi_check);
			projection_T0i_project(pcls_cdm, Bi_check, phi);
			if (sim.baryon_flag)
				projection_T0i_project(pcls_b, Bi_check, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_T0i_project(pcls_ncdm+i, Bi_check, phi);
			projection_T0i_comm(Bi_check);
			plan_Bi_check->execute(FFT_FORWARD);
			projectFTvector(*BiFT_check, *BiFT_check, fourpiG / (double) sim.numpts / (double) sim.numpts);
		}
		plan_Bi_check->execute(FFT_BACKWARD);
			
		for (x.first(); x.test(); x.next())
		{
			(*Bi_check)(x,0) /= a * a * sim.numpts;
			(*Bi_check)(x,1) /= a * a * sim.numpts;
			(*Bi_check)(x,2) /= a * a * sim.numpts;
		}
#ifdef EXTERNAL_IO
		Bi_check->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			Bi_check->saveHDF5_coarseGrain3D(h5filename + filename + "_B_check.h5", sim.downgrade_factor);
		else
			Bi_check->saveHDF5(h5filename + filename + "_B_check.h5");
#endif
	}
#endif

	if (sim.out_snapshot & MASK_GADGET)
	{
		hdr.time = a;
		hdr.redshift = (1./a) - 1.;
				
		hdr.npart[1] = (unsigned int) (sim.numpcl[0] / sim.tracer_factor[0]);
		if (sim.numpcl[0] % sim.tracer_factor[0]) hdr.npart[1]++;
		hdr.npartTotal[1] = hdr.npart[1];
		if (sim.baryon_flag)
			hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * cosmo.Omega_cdm * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
		else
			hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
		pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.tracer_factor[0]);
				
		if (sim.baryon_flag)
		{
			hdr.npart[1] = (unsigned int) (sim.numpcl[1] / sim.tracer_factor[1]);
			if (sim.numpcl[1] % sim.tracer_factor[1]) hdr.npart[1]++;
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
			pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.tracer_factor[1]);
		}
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			hdr.npart[1] = (unsigned int) (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]);
			if (sim.numpcl[i+1+sim.baryon_flag] % sim.tracer_factor[i+1+sim.baryon_flag]) hdr.npart[1]++;
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
			pcls_ncdm[i].saveGadget2(h5filename + filename + buffer, hdr, sim.tracer_factor[i+1+sim.baryon_flag]);
		}
	}
			
	if (sim.out_snapshot & MASK_PCLS)
	{
#ifdef EXTERNAL_IO
		pcls_cdm->saveHDF5_server_write();
		if (sim.baryon_flag)
			pcls_b->saveHDF5_server_write();
		for (i = 0; i < cosmo.num_ncdm; i++)
			pcls_ncdm[i].saveHDF5_server_write();
#else
		pcls_cdm->saveHDF5(h5filename + filename + "_cdm", 1);
		if (sim.baryon_flag)
			pcls_b->saveHDF5(h5filename + filename + "_b", 1);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			pcls_ncdm[i].saveHDF5(h5filename + filename + buffer, 1);
		}
#endif
	}
			
#ifdef EXTERNAL_IO
	ioserver.closeOstream();
#endif
}


//////////////////////////
// writeSpectra
//////////////////////////
// Description:
//   output of spectra
// 
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   pkcount        spectrum output index
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
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
//   Bi_check       pointer to allocated field (or NULL)
//   BiFT_check     pointer to allocated field (or NULL)
//   plan_Bi_check  pointer to FFT planner (or NULL)
//
// Returns:
// 
//////////////////////////

void writeSpectra(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const int pkcount, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, Field<Real> * Bi_check = NULL, Field<Cplx> * BiFT_check = NULL, PlanFFT<Cplx> * plan_Bi_check = NULL)
{
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[64];
	int i, j;
	Site x(phi->lattice());
	rKSite kFT(scalarFT->lattice());
    long numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;
	Cplx tempk;
	double Omega_ncdm;

	Real * kbin;
	Real * power;
	Real * kscatter;
	Real * pscatter;
	int * occupation;

	kbin = (Real *) malloc(sim.numbins * sizeof(Real));
	power = (Real *) malloc(sim.numbins * sizeof(Real));
	kscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	pscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	occupation = (int *) malloc(sim.numbins * sizeof(int));

	if (sim.out_pk & MASK_RBARE || sim.out_pk & MASK_DBARE || sim.out_pk & MASK_POT || ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag == 0))
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if (sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		for (i = 0; i < cosmo.num_ncdm; i++)
			scalarProjectionCIC_project(pcls_ncdm+i, source);
		scalarProjectionCIC_comm(source);
		plan_source->execute(FFT_FORWARD);
				
		if (sim.out_pk & MASK_RBARE || sim.out_pk & MASK_DBARE || ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag == 0))
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				
		if (sim.out_pk & MASK_RBARE)
		{
			sprintf(filename, "%s%s%03d_rhoN.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of rho_N", a);
		}

		if (sim.out_pk & MASK_DBARE)
		{
			sprintf(filename, "%s%s%03d_deltaN.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_m * cosmo.Omega_m, filename, "power spectrum of delta_N", a);
		}
				
		if (sim.out_pk & MASK_T00 && sim.gr_flag == 0)
		{
			sprintf(filename, "%s%s%03d_T00.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00", a);
		}

		if (sim.out_pk & MASK_DELTA && sim.gr_flag == 0)
		{
			sprintf(filename, "%s%s%03d_delta.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_m * cosmo.Omega_m, filename, "power spectrum of delta", a);
		}
				
		if (sim.out_pk & MASK_POT)
		{
			solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
			sprintf(filename, "%s%s%03d_psiN.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of psi_N", a);
		}
				
		if ((cosmo.num_ncdm > 0 || sim.baryon_flag) && (sim.out_pk & MASK_DBARE || (sim.out_pk & MASK_DELTA && sim.gr_flag == 0)))
		{
			projection_init(source);
			scalarProjectionCIC_project(pcls_cdm, source);
			scalarProjectionCIC_comm(source);
			plan_source->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			sprintf(filename, "%s%s%03d_cdm.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (sim.baryon_flag ? (cosmo.Omega_cdm * cosmo.Omega_cdm) : ((cosmo.Omega_cdm + cosmo.Omega_b) * (cosmo.Omega_cdm + cosmo.Omega_b))), filename, "power spectrum of delta_N for cdm", a);
			if (sim.baryon_flag)
			{
				// store k-space information for cross-spectra using SijFT as temporary array
				if (sim.out_pk & MASK_XSPEC)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, 0) = (*scalarFT)(kFT);
				}
				projection_init(source);
				scalarProjectionCIC_project(pcls_b, source);
				scalarProjectionCIC_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_b.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_b, filename, "power spectrum of delta_N for baryons", a);
				if (sim.out_pk & MASK_XSPEC)
				{
					extractCrossSpectrum(*scalarFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
					sprintf(filename, "%s%s%03d_cdmxb.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_cdm * cosmo.Omega_b, filename, "cross power spectrum of delta_N for cdm x baryons", a);
				}
			}
			Omega_ncdm = 0.;
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				projection_init(source);
				scalarProjectionCIC_project(pcls_ncdm+i, source);
				scalarProjectionCIC_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_ncdm%d.dat", sim.output_path, sim.basename_pk, pkcount, i);
				sprintf(buffer, "power spectrum of delta_N for ncdm %d", i);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_ncdm[i] * cosmo.Omega_ncdm[i], filename, buffer, a);
				Omega_ncdm += cosmo.Omega_ncdm[i];
				// store k-space information for cross-spectra using SijFT as temporary array
				if (cosmo.num_ncdm > 1 && i < 6)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, i) = (*scalarFT)(kFT);
				}						
			}
			if (cosmo.num_ncdm > 1 && cosmo.num_ncdm <= 7)
			{
				for (kFT.first(); kFT.test(); kFT.next())
				{
					for (i = 0; i < cosmo.num_ncdm-1; i++)
						(*scalarFT)(kFT) += (*SijFT)(kFT, i);
				}
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_ncdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * Omega_ncdm * Omega_ncdm, filename, "power spectrum of delta_N for total ncdm", a);
			}
			if (cosmo.num_ncdm > 1)
			{
				for (i = 0; i < cosmo.num_ncdm-1 && i < 5; i++)
				{
					for (j = i+1; j < cosmo.num_ncdm && j < 6; j++)
					{
						if (sim.out_pk & MASK_XSPEC || (i == 0 && j == 1) || (i == 2 && j == 3) || (i == 4 && j == 5))
						{
							extractCrossSpectrum(*SijFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR, i, j);
							sprintf(filename, "%s%s%03d_ncdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount, i, j);
							sprintf(buffer, "cross power spectrum of delta_N for ncdm %d x %d", i, j);
							writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_ncdm[i] * cosmo.Omega_ncdm[j], filename, buffer, a);
						}
					}
				}
			}
		}
	}
	
	if (sim.out_pk & MASK_PHI)
	{
		plan_phi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_phi.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi", a);
	}
			
	if (sim.out_pk & MASK_CHI)
	{
		plan_chi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_chi.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of chi", a);
	}
			
	if (sim.out_pk & MASK_HIJ)
	{
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if (sim.baryon_flag)
			projection_Tij_project(pcls_b, Sij, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_Tij_project(pcls_ncdm+i, Sij, a, phi);
		projection_Tij_comm(Sij);

		prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / (double) sim.numpts / (double) sim.numpts / a);
		plan_Sij->execute(FFT_FORWARD);
		projectFTtensor(*SijFT, *SijFT);

		extractPowerSpectrum(*SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_hij.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, 2. * M_PI * M_PI, filename, "power spectrum of hij", a);
	}
			
	if ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag > 0)
	{
		projection_init(source);
		projection_T00_project(pcls_cdm, source, a, phi);
		if (sim.baryon_flag)
			projection_T00_project(pcls_b, source, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_T00_project(pcls_ncdm+i, source, a, phi);
		projection_T00_comm(source);

		plan_source->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
		
		if (sim.out_pk & MASK_T00)
		{
			sprintf(filename, "%s%s%03d_T00.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00", a);
		}

		if (sim.out_pk & MASK_DELTA)
		{
			sprintf(filename, "%s%s%03d_delta.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)), filename, "power spectrum of delta", a);
		}
				
		if (cosmo.num_ncdm > 0 || sim.baryon_flag)
		{
			projection_init(source);
			projection_T00_project(pcls_cdm, source, a, phi);
			projection_T00_comm(source);
			plan_source->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			if (sim.out_pk & MASK_T00)
			{
				sprintf(filename, "%s%s%03d_T00cdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for cdm", a);
			}
			if (sim.out_pk & MASK_DELTA)
			{
				sprintf(filename, "%s%s%03d_deltacdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (sim.baryon_flag ? (cosmo.Omega_cdm * cosmo.Omega_cdm) : ((cosmo.Omega_cdm + cosmo.Omega_b) * (cosmo.Omega_cdm + cosmo.Omega_b))), filename, "power spectrum of delta for cdm", a);
			}
			if (sim.baryon_flag)
			{
				// store k-space information for cross-spectra using SijFT as temporary array
				if (sim.out_pk & MASK_XSPEC)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, 0) = (*scalarFT)(kFT);
				}
				projection_init(source);
				projection_T00_project(pcls_b, source, a, phi);
				projection_T00_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if (sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s%03d_T00b.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for baryons", a);
				}
				if (sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s%03d_deltab.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_b, filename, "power spectrum of delta for baryons", a);
				}
				if (sim.out_pk & MASK_XSPEC)
				{
					extractCrossSpectrum(*scalarFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
					sprintf(filename, "%s%s%03d_deltacdmxb.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_cdm, filename, "cross power spectrum of delta for cdm x baryons", a);
				}
			}
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				projection_init(source);
				projection_T00_project(pcls_ncdm+i, source, a, phi);
				projection_T00_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if (sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s%03d_T00ncdm%d.dat", sim.output_path, sim.basename_pk, pkcount, i);
					sprintf(buffer, "power spectrum of T00 for ncdm %d", i);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, buffer, a);
				}
				if (sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s%03d_deltancdm%d.dat", sim.output_path, sim.basename_pk, pkcount, i);
					sprintf(buffer, "power spectrum of delta for ncdm %d", i);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo, i) * bg_ncdm(a, cosmo, i), filename, buffer, a);
				}					
				// store k-space information for cross-spectra using SijFT as temporary array
				if (cosmo.num_ncdm > 1 && i < 6)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, i) = (*scalarFT)(kFT);
				}
			}
			if (cosmo.num_ncdm > 1 && cosmo.num_ncdm <= 7)
			{
				for (kFT.first(); kFT.test(); kFT.next())
				{
					for (i = 0; i < cosmo.num_ncdm-1; i++)
						(*scalarFT)(kFT) += (*SijFT)(kFT, i);
				}
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if (sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s%03d_T00ncdm.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for total ncdm", a);
				}
				if (sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s%03d_deltancdm.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo) * bg_ncdm(a, cosmo), filename, "power spectrum of delta for total ncdm", a);
				}
			}
			if (cosmo.num_ncdm > 1)
			{
				for (i = 0; i < cosmo.num_ncdm-1 && i < 5; i++)
				{
					for (j = i+1; j < cosmo.num_ncdm && j < 6; j++)
					{
						if (sim.out_pk & MASK_XSPEC || (i == 0 && j == 1) || (i == 2 && j == 3) || (i == 4 && j == 5))
						{
							extractCrossSpectrum(*SijFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR, i, j);
							if (sim.out_pk & MASK_T00)
							{
								sprintf(filename, "%s%s%03d_T00ncdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount, i, j);
								sprintf(buffer, "cross power spectrum of T00 for ncdm %d x %d", i, j);
								writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, buffer, a);
							}
							if (sim.out_pk & MASK_DELTA)
							{
								sprintf(filename, "%s%s%03d_deltancdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount, i, j);
								sprintf(buffer, "cross power spectrum of delta for ncdm %d x %d", i, j);
								writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo, i) * bg_ncdm(a, cosmo, j), filename, buffer, a);
							}
						}
					}
				}
			}
		}
	}
			
	if (sim.out_pk & MASK_B)
	{
		extractPowerSpectrum(*BiFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_B.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a);
			
#ifdef CHECK_B
		if (sim.vector_flag == VECTOR_PARABOLIC)
		{
			projection_init(Bi_check);
			projection_T0i_project(pcls_cdm, Bi_check, phi);
			if (sim.baryon_flag)
				projection_T0i_project(pcls_b, Bi_check, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_T0i_project(pcls_ncdm+i, Bi_check, phi);
			projection_T0i_comm(Bi_check);
			plan_Bi_check->execute(FFT_FORWARD);
			projectFTvector(*BiFT_check, *BiFT_check, fourpiG / (double) sim.numpts / (double) sim.numpts);
		}
		extractPowerSpectrum(*BiFT_check, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_B_check.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a);
#endif
	}

	free(kbin);
	free(power);
	free(kscatter);
	free(pscatter);
	free(occupation);
}

#endif

