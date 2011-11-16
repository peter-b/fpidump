/*
 * fpidump -- Extract data from Astrium .fpi files
 * Copyright (C) 2011  Peter Brett <p.brett@surrey.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include <tiffio.h>
#include <glib.h>

/* Representation of FPI file header */
typedef struct _FpiHeader FpiHeader;

struct _FpiHeader {
  char jobName[500];
  char rawDataFilename[500];
  char modeName[500];
  char priNamesProcessed[500];
  char subBandName[500];
  int txPolarisation;
  int channelProcessed;
  double imageStartPosition;
  double imageLength;
  char productName[500];
  int floatsPerPixel;
  int numberOfRows;
  int numberOfColumns;
  int slantRangeOrGroundRange;
  double azimuthSpatialResolution;
  double rangeSpatialResolution;
  double azimuthPixelSpacing;
  double rangePixelSpacing;
  double meanAltitude;
  double meanTerrainHeight;
  double minSlantRange;
  double numAzimuthLooks;
  double numRangeLooks;

  double firstPixelLongitude;
  double firstPixelLatitude;
  double lastPixelLongitude;
  double lastPixelLatitude;
  double imageOrientation;

  double acquisitionStartGPStimeOfWeek;
  double acquisitionEndGPStimeOfWeek;
  char absoluteDateOfAcquisition[500];

  bool spareBool[50];
  int spareInt[50];
  double spareDouble[50];
};

/* Position in FPI file of start of header */
#define FPI_HEADER_OFFSET 0

/* Position in FPI file of start of data */
#define FPI_DATA_OFFSET 4318

#define GETOPT_OPTIONS "m:t:h"

/* Help message */
static void
usage (char *name, int status)
{
  printf (
"Usage: %s OPTION... FPIFILE\n"
"\n"
"Extract SAR data from Astrium .fpi files.\n"
"\n"
"  -m METADATA     Export metadata to parameter file.\n"
"  -t TIFF         Export image to TIFF.\n"
"  -h              Display this message and exit.\n"
"\n"
"Extract image data and associated metadata from FPIFILE.  Metadata is\n"
"exported in an configuration file-like format, and image raster data\n"
"is exported as 32-bit floating point TIFF format.\n"
"\n"
"Please report bugs to %s.\n",
name, PACKAGE_BUGREPORT);
  exit (status);
}

/* Helper function for setting strings in a key file with a maximum
 * field size. */
static void
key_file_set_stringn (GKeyFile *key_file, const gchar *group_name,
                      const gchar *key, const gchar *str,
                      gsize len)
{
  gchar *buf = g_strndup (str, len);
  g_key_file_set_string (key_file, group_name, key, buf);
  g_free (buf);
}

/* Helper function that wraps fread(3) and quits with an error message
 * on failure.  The filename argument is used for generating error
 * messages. */
static void
fread_error (void *buf, size_t size, FILE *fp, const char *filename)
{
  size_t status = fread (buf, size, 1, fp);
  if (status != 1) {
    const char *msg;
    if (ferror (fp)) {
      msg = strerror (errno);
    } else if (feof (fp)) {
      msg = "Unexpected end-of-file";
    } else {
      msg = "Unexpected error";
    }
    fprintf (stderr, "ERROR: Could not read FPI header from %s: %s\n",
             filename, msg);
    exit (2);
  }
}

/* Read an FPI header from fpi_fpi.  The filename argument is used for
 * generating error messages. */
static void
read_header (FpiHeader *header, FILE *fpi_fp, const char *filename)
{
  fseek (fpi_fp, FPI_HEADER_OFFSET, SEEK_SET);

#define FREAD_ERROR(field, size) fread_error (&header->field, size, fpi_fp, filename)
  FREAD_ERROR (jobName, 500);
  FREAD_ERROR (rawDataFilename, 500);
  FREAD_ERROR (modeName, 500);
  FREAD_ERROR (priNamesProcessed, 500);
  FREAD_ERROR (subBandName, 500);
  FREAD_ERROR (txPolarisation, sizeof (int));
  FREAD_ERROR (channelProcessed, sizeof (int));
  FREAD_ERROR (imageStartPosition, sizeof (double));
  FREAD_ERROR (imageLength, sizeof (double));
  FREAD_ERROR (productName, 500);
  FREAD_ERROR (floatsPerPixel, sizeof (int));
  FREAD_ERROR (numberOfRows, sizeof (int));
  FREAD_ERROR (numberOfColumns, sizeof (int));
  FREAD_ERROR (slantRangeOrGroundRange, sizeof (int));
  FREAD_ERROR (azimuthSpatialResolution, sizeof (double));
  FREAD_ERROR (rangeSpatialResolution, sizeof (double));
  FREAD_ERROR (azimuthPixelSpacing, sizeof (double));
  FREAD_ERROR (rangePixelSpacing, sizeof (double));
  FREAD_ERROR (meanAltitude, sizeof (double));
  FREAD_ERROR (meanTerrainHeight, sizeof (double));
  FREAD_ERROR (minSlantRange, sizeof (double));
  FREAD_ERROR (numAzimuthLooks, sizeof (double));
  FREAD_ERROR (numRangeLooks, sizeof (double));

  FREAD_ERROR (firstPixelLongitude, sizeof (double));
  FREAD_ERROR (firstPixelLatitude, sizeof (double));
  FREAD_ERROR (lastPixelLongitude, sizeof (double));
  FREAD_ERROR (lastPixelLatitude, sizeof (double));
  FREAD_ERROR (imageOrientation, sizeof (double));

  FREAD_ERROR (acquisitionStartGPStimeOfWeek, sizeof (double));
  FREAD_ERROR (acquisitionEndGPStimeOfWeek, sizeof (double));

  FREAD_ERROR (spareBool, 50 * sizeof (bool));
  FREAD_ERROR (spareInt, 50 * sizeof (int));
  FREAD_ERROR (spareDouble, 50 * sizeof (double));
#undef FREAD_ERROR
}

/* Generate a metadata file for an FPI header, outputting to
 * filename. */
static void
generate_metadata_file (FpiHeader *header, const char *filename)
{
  GKeyFile *metadata;
  GError *err = NULL;

  /* Generate metadata keyfile */
  metadata = g_key_file_new ();

  key_file_set_stringn (metadata, "What", "JobName",
                        header->jobName, 500);
  g_key_file_set_comment (metadata, "What", "JobName",
                          "Name of the processing job", NULL);

  key_file_set_stringn (metadata, "What", "RawDataFilename",
                        header->rawDataFilename, 500);
  g_key_file_set_comment (metadata, "What", "RawDataFilename",
                          "Raw data filename", NULL);

  key_file_set_stringn (metadata, "What", "ModeName",
                        header->modeName, 500);
  g_key_file_set_comment (metadata, "What", "ModeName",
                          "Mode name", NULL);

  key_file_set_stringn (metadata, "What", "PRINamesProcessed",
                        header->priNamesProcessed, 500);
  g_key_file_set_comment (metadata, "What", "PRINamesProcessed",
                          "PRI names processed", NULL);

  key_file_set_stringn (metadata, "What", "SubBandName",
                        header->subBandName, 500);
  g_key_file_set_comment (metadata, "What", "SubBandName",
                          "Sub-band name", NULL);

  g_key_file_set_integer (metadata, "What", "TxPolarisation",
                          header->txPolarisation);
  g_key_file_set_comment (metadata, "What", "TxPolarisation",
                          "Transmit polarisation (V = 1, H = 2)", NULL);

  g_key_file_set_integer (metadata, "What", "ChannelProcessed",
                          header->channelProcessed);
  g_key_file_set_comment (metadata, "What", "ChannelProcessed",
                          "Channel processed, i.e. receive polarisation "
                          "(V = 1, H = 2)", NULL);

  g_key_file_set_double (metadata, "What", "ImageStartPosition",
                         header->imageStartPosition);
  g_key_file_set_comment (metadata, "What", "ImageStartPosition",
                          "Along-track position of image start relative to start "
                          "of whole image (m)", NULL);

  g_key_file_set_double (metadata, "What", "ImageLength",
                         header->imageLength);
  g_key_file_set_comment (metadata, "What", "ImageLength",
                          "Along-track length of image (m)", NULL);

  key_file_set_stringn (metadata, "What", "ProductName",
                        header->productName, 500);
  g_key_file_set_comment (metadata, "What", "ProductName",
                          "Name of product to be generated", NULL);

  g_key_file_set_integer (metadata, "What", "FloatsPerPixel",
                          header->floatsPerPixel);
  g_key_file_set_comment (metadata, "What", "FloatsPerPixel",
                          "32-bit floats per pixel", NULL);

  g_key_file_set_integer (metadata, "What", "NumberOfRows",
                          header->numberOfRows);
  g_key_file_set_comment (metadata, "What", "NumberOfRows",
                          "Number of pixel rows", NULL);

  g_key_file_set_integer (metadata, "What", "NumberOfColumns",
                          header->numberOfColumns);
  g_key_file_set_comment (metadata, "What", "NumberOfColumns",
                          "Number of pixel columns", NULL);

  g_key_file_set_integer (metadata, "What", "SlantRangeOrGroundRange",
                          header->slantRangeOrGroundRange);
  g_key_file_set_comment (metadata, "What", "SlantRangeOrGroundRange",
                          "Slant or ground range image (slant = 0, ground = 1)",
                          NULL);

  g_key_file_set_double (metadata, "What", "AzimuthSpatialResolution",
                         header->azimuthSpatialResolution);
  g_key_file_set_comment (metadata, "What", "AzimuthSpatialResolution",
                          "Azimuth spatial resolution (m)", NULL);

  g_key_file_set_double (metadata, "What", "RangeSpatialResolution",
                         header->rangeSpatialResolution);
  g_key_file_set_comment (metadata, "What", "RangeSpatialResolution",
                          "Range spatial resolution (m)", NULL);

  g_key_file_set_double (metadata, "What", "AzimuthPixelSpacing",
                         header->azimuthPixelSpacing);
  g_key_file_set_comment (metadata, "What", "AzimuthPixelSpacing",
                          "Azimuth pixel spacing (m)", NULL);

  g_key_file_set_double (metadata, "What", "RangePixelSpacing",
                         header->rangePixelSpacing);
  g_key_file_set_comment (metadata, "What", "RangePixelSpacing",
                          "Range pixel spacing (m)", NULL);

  g_key_file_set_double (metadata, "What", "MeanAltitude",
                         header->meanAltitude);
  g_key_file_set_comment (metadata, "What", "MeanAltitude",
                          "Mean altitude (m) above terrain assuming terrain is at "
                          "MeanTerrainHeight above WGS84 ellipsoid", NULL);

  g_key_file_set_double (metadata, "What", "MeanTerrainHeight",
                         header->meanTerrainHeight);
  g_key_file_set_comment (metadata, "What", "MeanTerrainHeight",
                          "Mean terrain height above WGS84 ellipsoid (m)", NULL);

  g_key_file_set_double (metadata, "What", "MinSlantRange",
                         header->minSlantRange);
  g_key_file_set_comment (metadata, "What", "MinSlantRange",
                          "Near-edge slant range (m)", NULL);

  g_key_file_set_double (metadata, "What", "NumAzimuthLooks",
                         header->numAzimuthLooks);
  g_key_file_set_comment (metadata, "What", "NumAzimuthLooks",
                          "Number of azimuth looks (1.0 for SLC images)", NULL);

  g_key_file_set_double (metadata, "What", "NumRangeLooks",
                         header->numRangeLooks);
  g_key_file_set_comment (metadata, "What", "NumRangeLooks",
                          "Number of range looks (1.0 for SLC images)", NULL);

  g_key_file_set_double (metadata, "Where", "FirstPixelLongitude",
                         header->firstPixelLongitude);
  g_key_file_set_comment (metadata, "Where", "FirstPixelLongitude",
                          "Approx. longitude of first pixel (degrees East)",
                          NULL);

  g_key_file_set_double (metadata, "Where", "FirstPixelLatitude",
                         header->firstPixelLatitude);
  g_key_file_set_comment (metadata, "Where", "FirstPixelLatitude",
                          "Approx. longitude of first pixel (degrees North)",
                          NULL);

  g_key_file_set_double (metadata, "Where", "LastPixelLongitude",
                         header->firstPixelLongitude);
  g_key_file_set_comment (metadata, "Where", "LastPixelLongitude",
                          "Approx. longitude of last pixel (degrees East)",
                          NULL);

  g_key_file_set_double (metadata, "Where", "LastPixelLatitude",
                         header->firstPixelLatitude);
  g_key_file_set_comment (metadata, "Where", "LastPixelLatitude",
                          "Approx. longitude of last pixel (degrees North)",
                          NULL);

  g_key_file_set_double (metadata, "Where", "ImageOrientation",
                         header->imageOrientation);
  g_key_file_set_comment (metadata, "Where", "ImageOrientation",
                          "Image orientation (degrees, North = 0)", NULL);

  g_key_file_set_double (metadata, "When", "AcquisitionStartGPSTimeOfWeek",
                         header->acquisitionStartGPStimeOfWeek);
  g_key_file_set_comment (metadata, "When", "AcquisitionStartGPSTimeOfWeek",
                          "Start time of entire acquisition since preceeding "
                          "Sunday 00:00", NULL);

  g_key_file_set_double (metadata, "When", "AcquisitionEndGPSTimeOfWeek",
                         header->acquisitionEndGPStimeOfWeek);
  g_key_file_set_comment (metadata, "When", "AcquisitionEndGPSTimeOfWeek",
                          "End time of entire acquisition since preceeding "
                          "Sunday 00:00", NULL);

  key_file_set_stringn (metadata, "When", "AbsoluteDateOfAcquisition",
                        header->absoluteDateOfAcquisition, 500);
  g_key_file_set_comment (metadata, "When", "AbsoluteDateOfAcquisition",
                          "Actual acquisition date", NULL);

  gchar *buf = g_key_file_to_data (metadata, NULL, NULL);
  gboolean status = g_file_set_contents (filename, buf, -1, &err);
  if (!status) {
    fprintf (stderr, "ERROR: Could not write metadata to %s: %s",
             filename, err->message);
    exit (3);
  }
  g_free (buf);
  g_key_file_free (metadata);
}

/* Generates a TIFF file named filename by reading from fpi_fp.  The
 * generated .TIFF file will be IEEEFP or COMPLEXIEEEFP depending on
 * the format of the source file. */
static void
generate_tiff_file (FpiHeader *header, FILE *fpi_fp, const char *filename)
{
  TIFF *tif = NULL;
  int rows_per_strip, num_strips;
  int result = 1;
  char *buffer = NULL;

  fseek (fpi_fp, FPI_DATA_OFFSET, SEEK_SET);

  /* Open TIFF file */
  tif = TIFFOpen (filename, "wb");
  if (!tif) {
    fprintf (stderr, "ERROR: Could not open TIFF file %s for writing.",
             filename);
    exit (4);
  }

  /* Set TIFF tags */
  TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, header->numberOfColumns);
  TIFFSetField (tif, TIFFTAG_IMAGELENGTH, header->numberOfRows);
  TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);

  if (header->floatsPerPixel == 2) {
    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 64);
    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT,
                  SAMPLEFORMAT_COMPLEXIEEEFP);

  } else if (header->floatsPerPixel == 1) {
    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT,
                  SAMPLEFORMAT_IEEEFP);
  }

  TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

  rows_per_strip = TIFFDefaultStripSize (tif, 0);
  TIFFSetField (tif, TIFFTAG_ROWSPERSTRIP, rows_per_strip);
  num_strips = TIFFNumberOfStrips (tif);

  /* Copy data from .FPI file into TIFF strips */
  size_t strip_size = TIFFStripSize (tif);
  size_t row_size = (header->numberOfColumns * header->floatsPerPixel
                     * sizeof (float));
  buffer = malloc (strip_size);

  for (int row = 0, strip = 0; strip < num_strips; strip++) {
    size_t len = 0;
    char *dest = buffer;

    /* Read rows from input file */
    for (int i = 0;
         (i < rows_per_strip) && (row < header->numberOfRows);
         i++, row++) {

      dest = buffer + i*row_size;
      result = fread (dest, row_size, 1, fpi_fp);
      len += row_size;

      if (result != 1) {
        const char *msg;
        if (ferror (fpi_fp)) {
          msg = strerror (errno);
        } else if (feof (fpi_fp)) {
          msg = "Unexpected end-of-file";
        } else {
          msg = "Unexpected error";
        }
        fprintf (stderr, "ERROR: Could not read image data: %s\n", msg);
        exit (2);
      }
    }

    /* Write rows to output file */
    result = (TIFFWriteEncodedStrip (tif, strip, buffer,
                                     len) != -1);
    if (!result) {
      fprintf (stderr, "ERROR: Could not write TIFF strip %i/%i to %s\n",
               strip+1, num_strips, filename);
      exit (4);
    }
  }

  /* Clean up */
  free (buffer);
  TIFFClose (tif);
}

int
main (int argc, char **argv)
{
  int c;
  char *infile = NULL;
  char *mdfile = NULL;
  char *tifffile = NULL;
  FpiHeader header;

  /* Parse command-line arguments */
  while ((c = getopt (argc, argv, GETOPT_OPTIONS)) != -1) {
    switch (c) {
    case 'm':
      mdfile = optarg;
      break;
    case 't':
      tifffile = optarg;
      break;
    case 'h':
      usage (argv[0], 0);
      break;
    case '?':
      if ((optopt != ':') && (strchr (GETOPT_OPTIONS, optopt) != NULL)) {
        fprintf (stderr, "ERROR: -%c option requires an argument.\n\n", optopt);
      } else if (isprint (optopt)) {
        fprintf (stderr, "ERROR: Unknown option -%c.\n\n", optopt);
      } else {
        fprintf (stderr, "ERROR: Unknown option character '\\x%x'.\n\n",
                 optopt);
      }
      usage (argv[0], 1);
      break;
    default:
      abort ();
    }
  }

  /* Get input and output filenames */
  if (argc - optind < 1) {
    fprintf (stderr, "ERROR: You must specify an input file.\n\n");
    usage (argv[0], 1);
  }
  infile = argv[optind];

  /* Attempt to load input file */
  FILE *fp = fopen (infile, "rb");
  if (fp == NULL) {
    const char *msg = errno ? strerror (errno) : "Unexpected error";
    fprintf (stderr, "ERROR: Could not load SAR data from %s: %s\n",
             infile, msg);
    exit (2);
  }

  /* Read metadata */
  read_header (&header, fp, infile);

  /* FIXME swab header */

  if (mdfile)
    generate_metadata_file (&header, mdfile);

  if (tifffile)
    generate_tiff_file (&header, fp, tifffile);

  fclose (fp);
  exit (0);
}
