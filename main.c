
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <getopt.h>

#include "nvutility.h"

#include "misp.h"
#include "chrtr2.h"

#include "version.h"


#define         FILTER 9
#define         EPS 1e-10


typedef struct
{
  CHRTR2_RECORD      ch2;
  int32_t            rank;
} CH2_GRID;



/*

    Programmer : Jan C. Depner
    Date : 01/18/11

    See usage (below) for an explanation of the purpose of this program.

*/


void usage ()
{
  fprintf (stderr, "\n\nUsage: chrtr2_merge [-e] [-b SIZE] [-n] CHRTR2_FILE1 CHRTR2_FILE2 [CHRTR2_FILE3...] [-o OUTPUT_FILE]\n\n");
  fprintf (stderr, "This program merges two or more CHRTR2 grids into a single CHRTR2 grid file.\n");
  fprintf (stderr, "The first file name on the command line takes precedence over the second\n");
  fprintf (stderr, "which takes precedence over the third... rinse, wash, repeat.  There is a\n");
  fprintf (stderr, "limit of 16 CHRTR2 files that can be merged.\n\n");
  fprintf (stderr, "-e = exclude\n");
  fprintf (stderr, "-b = buffer zone SIZE in grid cells for exclude (implies -e)\n");
  fprintf (stderr, "-n = no regrid of the output file\n");
  fprintf (stderr, "-o = set the output file name instead of defaulting\n\n");
  fprintf (stderr, "Examples:\n\n");
  fprintf (stderr, "chrtr2_merge file1.ch2 file2.ch2\n\n");
  fprintf (stderr, "  Inserts file1.ch2 into file1_merged.ch2.  Then inserts file2.ch2 into\n");
  fprintf (stderr, "  file1_merged.ch2 only where there is no data from file1.ch2.\n");
  fprintf (stderr, "  file1_merged.ch2 will be the same size and grid spacing as file1.ch2\n\n");
  fprintf (stderr, "chrtr2_merge file1.ch2 file2.ch2 file3.ch2\n\n");
  fprintf (stderr, "  This is the same as the first example except that it also inserts file3.ch2\n");
  fprintf (stderr, "  into file1_merged.ch2 where there is no data from file1.ch2 or file2.ch2\n\n");
  fprintf (stderr, "chrtr2_merge -n file1.ch2 file2.ch2\n\n");
  fprintf (stderr, "  This is the same as the first example except that file1_merged will not be\n");
  fprintf (stderr, "  regridded.\n\n");
  fprintf (stderr, "chrtr2_merge -n file1.ch2 file2.ch2 -o file1_file2_merged.ch2\n\n");
  fprintf (stderr, "  This is the same as the previous example except that the output file name\n");
  fprintf (stderr, "  will be file1_file2_merged.ch2 instead of file1_merged.ch2.\n\n");
  fprintf (stderr, "chrtr2_merge -e file1.ch2 file2.ch2\n\n");
  fprintf (stderr, "  Creates file1_merged.ch2 that has the grid spacing of file1.ch2 and has an\n");
  fprintf (stderr, "  MBR that includes both file1.ch2 and file2.ch2.  Data from file2.ch2 will be\n");
  fprintf (stderr, "  inserted anywhere there is no real data from file1.ch2 within four grid cells\n");
  fprintf (stderr, "  of the file2.ch2 data.\n\n");
  fprintf (stderr, "chrtr2_merge -b 10 file1.ch2 file2.ch2 file3.ch2\n\n");
  fprintf (stderr, "  Same as the above example except that the area of file1_merged.ch2 will be an\n");
  fprintf (stderr, "  MBR that includes all three files and the data from file3.ch2 will be\n");
  fprintf (stderr, "  inserted only where there are no points from file1.ch2 or file2.ch2 within 10\n");
  fprintf (stderr, "  grid cells of the data from file3.ch2.\n\n");

  fflush (stderr);
  exit (-1);
}



int32_t main (int32_t argc, char *argv[])
{
  char               c;
  extern char        *optarg;
  extern int         optind;
  int32_t            i, j, k, m, n, option_index = 0, chrtr2_handle[17], buffer_size = 4, start_x, end_x, start_y, end_y;
  int32_t            row_filter, col_filter, percent = 0 , old_percent = -1, grid_rows, grid_cols, input_count = 0, file_count = 0;
  char               input_file[16][512], output_file[512];
  uint8_t            exclude = NVFalse, dateline = NVFalse, regrid = NVTrue, hit;
  CHRTR2_HEADER      chrtr2_header[17];
  CHRTR2_RECORD      chrtr2_record;
  CH2_GRID           **grid;
  float              min_z, max_z, *array;
  double             lat, lon;
  NV_F64_MBR         new_mbr;
  NV_F64_XYMBR       mbr, misp_mbr;
  NV_F64_COORD3      xyz;
  NV_F64_COORD2      xy;
  NV_I32_COORD2      coord, coord2;


  printf ("\n\n %s \n\n\n", VERSION);


  strcpy (output_file, "");

  while (NVTrue) 
    {
      static struct option long_options[] = {{0, no_argument, 0, 0}};

      c = (char) getopt_long (argc, argv, "enb:o:", long_options, &option_index);
      if (c == -1) break;

      switch (c) 
        {
        case 0:

          switch (option_index)
            {
            case 0:
              break;
            }
          break;

        case 'e':
          exclude = NVTrue;
          break;

        case 'n':
          regrid = NVFalse;
          break;

        case 'b':
          sscanf (optarg, "%d", &buffer_size);
          exclude = NVTrue;
          break;

        case 'o':
          strcpy (output_file, optarg);
          break;

        default:
          usage ();
          exit (-1);
          break;
        }
    }


  /* Make sure we got the mandatory file names.  */

  file_count = argc - optind;
  if (file_count < 2 || file_count > 16) usage ();


  /*  Open all of the input files and determine the MBR of the output file.  */

  new_mbr.wlon = 999.0;
  new_mbr.elon = -999.0;
  new_mbr.slat = 999.0;
  new_mbr.nlat = -999.0;

  for (i = 0 ; i < file_count ; i++)
    {
      strcpy (input_file[i], argv[optind + i]);

      fprintf (stderr, "Input file %d  : %s\n", i + 1, input_file[i]);
      fflush (stderr);


      /*  Open the input file.  */

      chrtr2_handle[i] = chrtr2_open_file (input_file[i], &chrtr2_header[i], CHRTR2_READONLY);

      if (chrtr2_handle[i] < 0)
        {
          fprintf (stderr, "\n\nThe file %s is not a CHRTR2 file or there was an error reading the file.\nThe error message returned was:%s\n\n",
                   input_file[i], chrtr2_strerror ());
          exit (-1);
        }

      new_mbr.wlon = MIN (new_mbr.wlon, chrtr2_header[i].mbr.wlon);
      new_mbr.slat = MIN (new_mbr.slat, chrtr2_header[i].mbr.slat);
      new_mbr.elon = MAX (new_mbr.elon, chrtr2_header[i].mbr.elon);
      new_mbr.nlat = MAX (new_mbr.nlat, chrtr2_header[i].mbr.nlat);

      if (!dateline && new_mbr.elon > 360.0) dateline = NVTrue;
    }


  if (dateline && new_mbr.elon < new_mbr.wlon) new_mbr.elon += 360.0;


  chrtr2_header[16] = chrtr2_header[0];
  chrtr2_header[16].mbr = new_mbr;
  chrtr2_header[16].width = NINT ((new_mbr.elon - new_mbr.wlon) / chrtr2_header[0].lon_grid_size_degrees) + 1;
  chrtr2_header[16].height = NINT ((new_mbr.nlat - new_mbr.slat) / chrtr2_header[0].lat_grid_size_degrees) + 1;


  /*  Make the output file name.  */

  if (strlen (output_file) < 3)
    {
      strcpy (output_file, input_file[0]);
      sprintf (&output_file[strlen (output_file) - 4], "__merged.ch2");
    }
  else
    {
      /*  Make sure the .ch2 extension was included if the output file was specified on the command line.  */

      if (strcmp (&output_file[strlen (output_file) - 4], ".ch2")) strcat (output_file, ".ch2");
    }


  /*  Try to create and open the chrtr2 output file.  */

  chrtr2_handle[16] = chrtr2_create_file (output_file, &chrtr2_header[16]);
  if (chrtr2_handle[16] < 0)
    {
      chrtr2_perror ();
      exit (-1);
    }


  fprintf (stderr, "Output file : %s\n\n", output_file);
  fflush (stderr);


  /*  Allocate the output grid in memory so we don't have to keep reading and writing the output file.  */

  grid = (CH2_GRID **) malloc (chrtr2_header[16].height * sizeof (CH2_GRID *));
  if (grid == NULL)
    {
      perror ("Allocating grid array in main.c");
      exit (-1);
    }

  for (i = 0 ; i < chrtr2_header[16].height ; i++)
    {
      grid[i] = (CH2_GRID *) calloc (chrtr2_header[16].width, sizeof (CH2_GRID));
      if (grid[i] == NULL)
        {
          perror ("Allocating grid[i] array in main.c");
          exit (-1);
        }
    }
      

  /*  Read all of the input CHRTR2 files and fill the sparse grid.  */

  for (i = 0 ; i < file_count ; i++)
    {
      /*  Loop for height of input file.  */

      for (j = 0 ; j < chrtr2_header[i].height ; j++)
        {
          coord.y = j;


          /*  Loop for width of input file.  */

          for (k = 0 ; k < chrtr2_header[i].width ; k++)
            {
              coord.x = k;


              /*  Read the input record.  */

              chrtr2_read_record (chrtr2_handle[i], coord, &chrtr2_record);


              /*  For the first file we just slap the data into the grid.  */

	      if (!i)
		{
		  grid[coord2.y][coord2.x].ch2 = chrtr2_record;
		  grid[coord2.y][coord2.x].rank = i + 1;
		}


              /*  Get the lat and lon of the center position of the input grid cell.  */

              chrtr2_get_lat_lon (chrtr2_handle[i], &lat, &lon, coord);

              lat = lat + EPS;
              lon = lon + EPS;      


              /*  Check for dateline crossing.  */

              if (dateline && lon < 0.0) lon += 360.0;


              /*  Check to see if the lat and lon is in the output file (it damn well should be).  */

              if (!chrtr2_get_coord (chrtr2_handle[16], lat, lon, &coord2))
                {
                  /*  If we're using the exclude option...  */

                  if (exclude)
                    {
                      /*  Check for real, hand-drawn/digitized, or land masked data.  */

                      if (chrtr2_record.status & (CHRTR2_REAL | CHRTR2_DIGITIZED_CONTOUR | CHRTR2_LAND_MASK))
                        {
                          /*  The first file is already populated (see above).  For subsequent files we have to check against what's already loaded into the grid.
                              Don't forget the buffer-size.  */

                          if (i)
                            {
                              /*  Determine the buffer_size box to exclude.  */

                              start_x = MAX (coord2.x - buffer_size, 0);
                              end_x = MIN (coord2.x + buffer_size, chrtr2_header[16].width - 1);
                              start_y = MAX (coord2.y - buffer_size, 0);
                              end_y = MIN (coord2.y + buffer_size, chrtr2_header[16].height - 1);


                              /*  Check all bins in the buffer.  */

                              hit = NVFalse;
                              for (m = start_y ; m <= end_y ; m++)
                                {
                                  for (n = start_x ; n <= end_x ; n++)
                                    {
                                      /*  First check to see that the point isn't from this file.  */

                                      if ((grid[m][n].rank != i + 1) && 
                                          (grid[m][n].ch2.status & (CHRTR2_REAL | CHRTR2_DIGITIZED_CONTOUR | CHRTR2_LAND_MASK))) hit = NVTrue;
                                    }
                                }


                              /*  If no bins in the buffer had real, hand-drawn/digitized, or land mask data, fill the bin.  */

                              if (!hit)
                                {
                                  grid[coord2.y][coord2.x].ch2 = chrtr2_record;
                                  grid[coord2.y][coord2.x].rank = i + 1;
                                }                                  
                            }
                        }
                    }
                  else
                    {
                      /*  The first file is already populated (see above).  For subsequent files we have to check against what's already loaded into the grid.
                          Don't forget the buffer-size.  */

                      if (i)
                        {
                          /*  We only load data where there is no data (i.e. NULL).  This is actually more of an insert than a merge but
                              this is what we need.  If there is no data in the bin, place the new data into the bin.  */

                          if (!grid[coord2.y][coord2.x].ch2.status)
                            {
                              grid[coord2.y][coord2.x].ch2 = chrtr2_record;
                              grid[coord2.y][coord2.x].rank = i + 1;
                            }
                        }
                    }
                }
            }


          percent = NINT (((float) j / (float) chrtr2_header[i].height) * 100.0);
          if (percent != old_percent)
            {
              fprintf (stderr, "Reading CHRTR2 file %d of %d - %03d%% complete\r", i + 1, file_count, percent);
              fflush (stderr);
              old_percent = percent;
            }
        }
    }

  fprintf (stderr, "                                                                   \r");
  fprintf (stderr, "\nData read complete\n\n");
  fflush (stderr);


  min_z = 9999999999.0;
  max_z = -9999999999.0;


  /*  Close the input files.  */

  for (i = 0 ; i < file_count ; i++) chrtr2_close_file (chrtr2_handle[i]);


  /*  Check to see if we want to regrid.  */

  if (regrid)
    {
      /*  Define the MBR for the new grid (adding the filter border).  */

      mbr.min_x = chrtr2_header[16].mbr.wlon;
      mbr.min_y = chrtr2_header[16].mbr.slat;
      mbr.max_x = chrtr2_header[16].mbr.elon;
      mbr.max_y = chrtr2_header[16].mbr.nlat;


      /*  Add the filter border to the MBR.  */

      mbr.min_x -= ((double) FILTER * chrtr2_header[16].lon_grid_size_degrees);
      mbr.min_y -= ((double) FILTER * chrtr2_header[16].lat_grid_size_degrees);
      mbr.max_x += ((double) FILTER * chrtr2_header[16].lon_grid_size_degrees);
      mbr.max_y += ((double) FILTER * chrtr2_header[16].lat_grid_size_degrees);


      /*  Number of rows and columns in the area  */

      grid_rows = NINT ((mbr.max_y - mbr.min_y) / chrtr2_header[16].lat_grid_size_degrees);
      grid_cols = NINT ((mbr.max_x - mbr.min_x) / chrtr2_header[16].lon_grid_size_degrees);


      row_filter = grid_rows - FILTER;
      col_filter = grid_cols - FILTER;


      /*  We're going to let MISP/SURF handle everything in zero based units of the bin size.  That is, we subtract off the
          west lon from longitudes then divide by the grid size in the X direction.  We do the same with the latitude using
          the south latitude.  This will give us values that range from 0.0 to grid5_cols in longitude and 0.0 to
          grid5_rows in latitude.  */

      misp_mbr.min_x = 0.0;
      misp_mbr.min_y = 0.0;
      misp_mbr.max_x = (double) grid_cols;
      misp_mbr.max_y = (double) grid_rows;


      misp_init (1.0, 1.0, 0.05, 4, 20.0, 20, 999999.0, -999999.0, -2, misp_mbr);


      for (i = 0 ; i < chrtr2_header[16].height ; i++)
        {
          coord.y = i;

          for (j = 0 ; j < chrtr2_header[16].width ; j++)
            {
              coord.x = j;


              /*  No point in loading null values.  */

              if (grid[coord.y][coord.x].ch2.status)
                {
                  chrtr2_get_lat_lon (chrtr2_handle[16], &xy.y, &xy.x, coord);


                  /*
                    Load the points.

                    IMPORTANT NOTE:  MISP and GMT (by default) grid using corner posts.  That is, the data in a bin is assigned to the 
                    lower left corner of the bin.  Normal gridding/binning systems use the center of the bin.  Because of this we need
                    to lie to MISP/GMT and tell them that the point is really half a bin lower and to the left.  This is extremely
                    confusing but it works ;-)
                  */

                  xyz.x = (xy.x - mbr.min_x) / chrtr2_header[16].lon_grid_size_degrees;
                  xyz.y = (xy.y - mbr.min_y) / chrtr2_header[16].lat_grid_size_degrees;
                  xyz.z = grid[coord.y][coord.x].ch2.z;

                  input_count++;

                  misp_load (xyz);
                }
            }

          percent = NINT (((float) i / (float) chrtr2_header[16].height) * 100.0);
          if (percent != old_percent)
            {
              fprintf (stderr, "Loading data for re-grid - %03d%% complete\r", percent);
              fflush (stderr);
              old_percent = percent;
            }
        }

      fprintf (stderr, "                                                                   \r");
      fprintf (stderr, "\nData load complete, %d points loaded\n\n", input_count);


      fprintf (stderr, "Processing grid\n");
      fflush (stderr);


      misp_proc ();


      fprintf (stderr, "Processing grid complete\n");
      fflush (stderr);


      array = (float *) malloc ((grid_cols + 1) * sizeof (float));

      if (array == NULL)
        {
          perror ("Allocating array in main.c");
          exit (-1);
        }


      /*  This is where we stuff the new interpolated surface into the new CHRTR2.  */

      for (i = 0 ; i < grid_rows ; i++)
        {
          if (!misp_rtrv (array)) break;


          /*  Only use data that aren't in the filter border  */

          if (i >= FILTER && i < row_filter)
            {
              coord.y = i - FILTER;

              for (j = 0 ; j < grid_cols ; j++)
                {
                  /*  Only use data that aren't in the filter border  */

                  if (j >= FILTER && j < col_filter)
                    {
                      coord.x = j - FILTER;


                      /*  Make sure we're inside the CHRTR2 bounds.  */

                      if (coord.y >= 0 && coord.y < chrtr2_header[16].height && coord.x >= 0 && coord.x < chrtr2_header[16].width)
                        {
                          /*  Don't replace real, hand-drawn/digitized, or land masked data.  */

                          if (!(grid[coord.y][coord.x].ch2.status & (CHRTR2_REAL | CHRTR2_DIGITIZED_CONTOUR | CHRTR2_LAND_MASK)))
                            {
                              grid[coord.y][coord.x].ch2.z = array[j];
                              grid[coord.y][coord.x].ch2.status |= CHRTR2_INTERPOLATED;
                            }

                          min_z = MIN (grid[coord.y][coord.x].ch2.z, min_z);
                          max_z = MAX (grid[coord.y][coord.x].ch2.z, max_z);

                          chrtr2_write_record (chrtr2_handle[16], coord, grid[coord.y][coord.x].ch2);
                        }
                    }
                }
            }

          percent = NINT (((float) i / (float) grid_rows) * 100.0);
          if (percent != old_percent)
            {
              fprintf (stderr, "Retrieving data for output file - %03d%% complete\r", percent);
              fflush (stderr);
              old_percent = percent;
            }
        }


      fprintf (stderr, "                                                                   \r");
      fprintf (stderr, "\nFinal grid retrieval complete\n\n");
      fflush (stderr);


      free (array);
    }
  else
    {
      for (i = 0 ; i < chrtr2_header[16].height ; i++)
        {
          coord.y = i;

          for (j = 0 ; j < chrtr2_header[16].width ; j++)
            {
              coord.x = j;

              if (grid[coord.y][coord.x].ch2.status)
                {
                  min_z = MIN (grid[coord.y][coord.x].ch2.z, min_z);
                  max_z = MAX (grid[coord.y][coord.x].ch2.z, max_z);

                  chrtr2_write_record (chrtr2_handle[16], coord, grid[coord.y][coord.x].ch2);
                }
            }

          percent = NINT (((float) i / (float) chrtr2_header[16].height) * 100.0);
          if (percent != old_percent)
            {
              fprintf (stderr, "Writing chrtr2 data - %03d%% complete\r", percent);
              fflush (stderr);
              old_percent = percent;
            }
        }


      fprintf (stderr, "                                                                   \r");
      fprintf (stderr, "\nFile writing complete\n\n");
      fflush (stderr);
    }


  for (i = 0 ; i < chrtr2_header[16].height ; i++) free (grid[i]);
  free (grid);

      
  chrtr2_close_file (chrtr2_handle[16]);


  /*  Update the header with the observed min and max values.  */

  chrtr2_header[16].min_observed_z = min_z;
  chrtr2_header[16].max_observed_z = max_z;

  chrtr2_update_header (chrtr2_handle[16], chrtr2_header[16]);

  chrtr2_close_file (chrtr2_handle[16]);


  fprintf (stderr, "\n\n%s complete\n\n\n", argv[0]);
  fflush (stderr);


  /*  Please ignore the following line.  It is useless.  Except...

      On some versions of Ubuntu, if I compile a program that doesn't use the math
      library but it calls a shared library that does use the math library I get undefined
      references to some of the math library functions even though I have -lm as the last
      library listed on the link line.  This happens whether I use qmake to build the
      Makefile or I have a pre-built Makefile.  Including math.h doesn't fix it either.
      The following line forces the linker to bring in the math library.  If there is a
      better solution please let me know at area.based.editor AT gmail DOT com.  */

  float ubuntu; ubuntu = 4.50 ; ubuntu = fmod (ubuntu, 1.0);


  return (0);
}
