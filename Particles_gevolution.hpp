//////////////////////////
// Particles_gevolution.hpp
//////////////////////////
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: November 2016
//
//////////////////////////

#ifndef PARTICLES_GEVOLUTION_HEADER
#define PARTICLES_GEVOLUTION_HEADER

#ifndef PCLBUFFER
#define PCLBUFFER 1048576
#endif

using namespace LATfield2;

template <typename part, typename part_info, typename part_dataType>
class Particles_gevolution: public Particles<part, part_info, part_dataType>
{
	public:
		void saveGadget2(string filename, gadget2_header & hdr, const int tracer_factor = 1);
		void loadGadget2(string filename, gadget2_header & hdr);
};

template <typename part, typename part_info, typename part_dataType>
void Particles_gevolution<part,part_info,part_dataType>::saveGadget2(string filename, gadget2_header & hdr, const int tracer_factor)
{
	float * posdata;
	float * veldata;
	void * IDs;
	MPI_File outfile;
	long count, npart;
	MPI_Offset offset_pos, offset_vel, offset_ID;
	MPI_Status status;
	uint32_t blocksize;
	uint32_t i;
	char fname[filename.length()+1];
	double rescale_vel = 1. / sqrt(hdr.time) / GADGET_VELOCITY_CONVERSION;
	
	filename.copy(fname, filename.length());
	fname[filename.length()] = '\0';
	
	LATfield2::Site xPart(this->lat_part_);
	typename std::list<part>::iterator it;
	
	if (hdr.num_files != 1)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": writing multiple Gadget2 files not currently supported!" << endl;
		return;
	}
	
	posdata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);
	veldata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);

#if GADGET_ID_BYTES == 8
	IDs = malloc(sizeof(int64_t) * PCLBUFFER);
#else
	IDs = malloc(sizeof(int32_t) * PCLBUFFER);
#endif
	
	npart = 0;
	for(xPart.first(); xPart.test(); xPart.next())
	{
		if(this->field_part_(xPart).size!=0)
		{
			for (it=(this->field_part_)(xPart).parts.begin(); it != (this->field_part_)(xPart).parts.end(); ++it)
			{
				if ((*it).ID % tracer_factor == 0)
					npart++;
			}
		}
	}
	
	if (parallel.rank() == 0)
	{
		parallel.send<long>(npart, 1);
		parallel.receive<long>(count, parallel.size()-1);
		if (count != hdr.npart[1]) cout << " error: number of particles in saveGadget2 does not match request!" << endl;
		count = 0;
	}
	else
	{
		parallel.receive<long>(count, parallel.rank()-1);
		npart += count;
		parallel.send<long>(npart, (parallel.rank()+1)%parallel.size());
	}
	
	MPI_File_open(parallel.lat_world_comm(), fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &outfile);
	
	offset_pos = (MPI_Offset) hdr.npart[1];
	offset_pos *= (MPI_Offset) (6 * sizeof(float) + ((GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t)));
	offset_pos += (MPI_Offset) (8 * sizeof(uint32_t) + sizeof(hdr));
	MPI_File_set_size(outfile, offset_pos);
	
	offset_pos = (MPI_Offset) (3 * sizeof(uint32_t) + sizeof(hdr)) + ((MPI_Offset) count) * ((MPI_Offset) (3 * sizeof(float)));
	offset_vel = offset_pos + (MPI_Offset) (2 * sizeof(uint32_t)) + ((MPI_Offset) hdr.npart[1]) * ((MPI_Offset) (3 * sizeof(float)));
	offset_ID = offset_vel + (MPI_Offset) (2 * sizeof(uint32_t)) + ((MPI_Offset) hdr.npart[1] - (MPI_Offset) count) * ((MPI_Offset) (3 * sizeof(float))) + ((MPI_Offset) count) * ((MPI_Offset) ((GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t)));
	
	if (parallel.rank() == 0)
	{
		blocksize = sizeof(hdr);		
		MPI_File_write_at(outfile, 0, &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, sizeof(uint32_t), &hdr, sizeof(hdr), MPI_BYTE, &status);
		MPI_File_write_at(outfile, sizeof(hdr) + sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		blocksize = 3 * sizeof(float) * hdr.npart[1];
		MPI_File_write_at(outfile, sizeof(hdr) + 2*sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_vel - 2*sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_vel - sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_ID - 2*sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		blocksize = ((GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t)) * hdr.npart[1];
		MPI_File_write_at(outfile, offset_ID - sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_ID + blocksize, &blocksize, 1, MPI_UNSIGNED, &status);
	}
	
	count = 0;
	for(xPart.first(); xPart.test(); xPart.next())
	{
		if(this->field_part_(xPart).size!=0)
		{
			for (it=(this->field_part_)(xPart).parts.begin(); it != (this->field_part_)(xPart).parts.end(); ++it)
			{
				if ((*it).ID % tracer_factor == 0)
				{
					for (i = 0; i < 3; i++)
						posdata[3*count+i] = (*it).pos[i] * hdr.BoxSize;
					
					for (i = 0; i < 3; i++)
						veldata[3*count+i] = (*it).vel[i] * rescale_vel / hdr.time;
					
#if GADGET_ID_BYTES == 8
					*((int64_t *) IDs + count) = (int64_t) (*it).ID;
#else	
					*((int32_t *) IDs + count) = (int32_t) (*it).ID;
#endif
					
					count++;
						
					if (count == PCLBUFFER)
					{
						MPI_File_write_at(outfile, offset_pos, posdata, 3 * count, MPI_FLOAT, &status);
						offset_pos += 3 * PCLBUFFER * sizeof(float);
						MPI_File_write_at(outfile, offset_vel, veldata, 3 * count, MPI_FLOAT, &status);
						offset_vel += 3 * PCLBUFFER * sizeof(float);
						count *= (GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t);
						MPI_File_write_at(outfile, offset_ID, IDs, count, MPI_BYTE, &status);
						offset_ID += count;
						count = 0;
					}
				}
			}
		}
	}
		
	if (count > 0)
	{
			MPI_File_write_at(outfile, offset_pos, posdata, 3 * count, MPI_FLOAT, &status);
			MPI_File_write_at(outfile, offset_vel, veldata, 3 * count, MPI_FLOAT, &status);
			count *= (GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t);
			MPI_File_write_at(outfile, offset_ID, IDs, count, MPI_BYTE, &status);
	}
	
	MPI_File_close(&outfile);
	
	free(posdata);
	free(veldata);
	free(IDs);
}


template <typename part, typename part_info, typename part_dataType>
void Particles_gevolution<part,part_info,part_dataType>::loadGadget2(string filename, gadget2_header & hdr)
{
	float * posdata;
	float * veldata;
	void * IDs;
	part pcl;
	MPI_File infile;
	uint32_t i, count, npart = 0;
	MPI_Offset offset_pos, offset_vel, offset_ID;
	MPI_Status status;
	uint32_t blocksize;
	char fname[filename.length()+1];
	double rescale_vel = 1. / GADGET_VELOCITY_CONVERSION;
	
	filename.copy(fname, filename.length());
	fname[filename.length()] = '\0';
	
	posdata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);
	veldata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);

#if GADGET_ID_BYTES == 8
	IDs = malloc(sizeof(int64_t) * PCLBUFFER);
#else
	IDs = malloc(sizeof(int32_t) * PCLBUFFER);
#endif

	MPI_File_open(parallel.lat_world_comm(), fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);

	MPI_File_read_all(infile, &blocksize, 1, MPI_UNSIGNED, &status);

	if (blocksize != sizeof(hdr))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": file type not recognized when reading Gadget2 file!" << endl;
		return;
	}

	MPI_File_read_all(infile, &hdr, sizeof(hdr), MPI_BYTE, &status);

	rescale_vel /= sqrt(hdr.time);
	offset_pos = (MPI_Offset) sizeof(hdr) + (MPI_Offset) (3 * sizeof(uint32_t));
	offset_vel = offset_pos + ((MPI_Offset) hdr.npart[1]) * ((MPI_Offset) (3 * sizeof(float))) + (MPI_Offset) (2 * sizeof(uint32_t));
	offset_ID = offset_vel + offset_vel - offset_pos;

	MPI_File_seek(infile, offset_pos, MPI_SEEK_SET);
	while (npart < hdr.npart[1])
	{
		count = (hdr.npart[1] - npart > PCLBUFFER) ? PCLBUFFER : (hdr.npart[1] - npart);

		MPI_File_read_all(infile, posdata, 3 * count, MPI_FLOAT, &status);
		offset_pos += (MPI_Offset) (3 * count * sizeof(float));
		MPI_File_seek(infile, offset_vel, MPI_SEEK_SET);
		MPI_File_read_all(infile, veldata, 3 * count, MPI_FLOAT, &status);
		offset_vel += (MPI_Offset) (3 * count * sizeof(float));
		MPI_File_seek(infile, offset_ID, MPI_SEEK_SET);
#if GADGET_ID_BYTES == 8
		MPI_File_read_all(infile, IDs, count * sizeof(int64_t), MPI_BYTE, &status);
		offset_ID += (MPI_Offset) (count * sizeof(int64_t));
#else
		MPI_File_read_all(infile, IDs, count * sizeof(int32_t), MPI_BYTE, &status);
		offset_ID += (MPI_Offset) (count * sizeof(int32_t));
#endif
		MPI_File_seek(infile, offset_pos, MPI_SEEK_SET);

		for (i = 0; i < 3 * count; i++)
		{
			posdata[i] /= hdr.BoxSize;
			veldata[i] *= hdr.time / rescale_vel;
		}

		for (i = 0; i < count; i++)
		{
#if GADGET_ID_BYTES == 8
			pcl.ID = *((int64_t *) IDs + i);
#else
			pcl.ID = *((int32_t *) IDs + i);
#endif
			pcl.pos[0] = posdata[3*i];
			pcl.pos[1] = posdata[3*i+1];
			pcl.pos[2] = posdata[3*i+2];
			pcl.vel[0] = veldata[3*i];
			pcl.vel[1] = veldata[3*i+1];
			pcl.vel[2] = veldata[3*i+2];
			this->addParticle_global(pcl);
		}
		
		npart += count;
	}

	MPI_File_close(&infile);
	
	free(posdata);
	free(veldata);
	free(IDs);
}

#endif
