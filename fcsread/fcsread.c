/*
 *  fcsread.c
 *  mpiStochEst
 *
 *	FCS file reader
 *
 *  This file is part of mpiStochEst.
 *  Copyright 2011-2017 Gabriele Lillacci.
 *
 *  mpiStochEst is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  mpiStochEst is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with mpiStochEst.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libxml2/libxml/tree.h>
#include <libxml2/libxml/parser.h>
#include <libxml2/libxml/xmlversion.h>
#include <libxml2/libxml/xmlstring.h>
#include <libxml2/libxml/xpath.h>

#include <parestlib.h>


int processPriorData (xmlDocPtr doc, double lb[], double ub[]);
int processDataSet (xmlDocPtr doc, gsl_vector * times, gsl_matrix * data);

// Global variable declarations
size_t M;	// Number of experimental samples to use for identification
size_t P;	// Number of measured species
size_t K;	// Number of time points
size_t R;	// Number of parameters


int main (int argc, char * argv[])
{
	// Initialize libxml2
	LIBXML_TEST_VERSION

	// Parse the problem file
	xmlDocPtr problem = xmlReadFile ("./fcsread/problem.xml", NULL, 0);
	if (problem == NULL)
	{
		fprintf (stderr, "Error reading problem file");
		return EXIT_FAILURE;
	}

	// Initialize XPath evaluation context
	xmlXPathContextPtr ctxt = xmlXPathNewContext (problem);

	// Initialize a xmlXPathObjectPtr object
	xmlXPathObjectPtr result;

	// Use XPath to find the number of samples to use for identification
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/samples", ctxt);
	M = (size_t) xmlXPathCastToNumber (result);
	xmlXPathFreeObject (result);

	// Use XPath to find the number of measured species
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/outputs", ctxt);
	P = (size_t) xmlXPathCastToNumber (result);
	xmlXPathFreeObject (result);

	// Use XPath to find the number of time points
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/timepoints", ctxt);
	K = (size_t) xmlXPathCastToNumber (result);
	xmlXPathFreeObject (result);

	// Use XPath to find the number of parameters
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/parameters", ctxt);
	R = (size_t) xmlXPathCastToNumber (result);
	xmlXPathFreeObject (result);

	// Use the results to initialize some useful variables
	double lb[R], ub[R];
	gsl_matrix * data = gsl_matrix_alloc (M, K*P);
	gsl_vector * times = gsl_vector_alloc (K);

	// Process the prior information
	if (processPriorData (problem, lb, ub) != GSL_SUCCESS)
	{
		fprintf (stderr, "error: the problem definition file does not contain valid prior information\n");
		return EXIT_FAILURE;
	}

	// Load the data from the FCS files
	if (processDataSet (problem, times, data) != GSL_SUCCESS)
	{
		fprintf (stderr, "error: could not process data set information\n");
		return EXIT_FAILURE;
	}

	print_gsl_vector (times);

	gsl_matrix_free (data);
	gsl_vector_free (times);
	printf ("\n");

	// Free XML resources
	xmlXPathFreeContext (ctxt);
	xmlFreeDoc (problem);

	// Shut down libxml2
	xmlCleanupParser ();
	xmlMemoryDump ();

	// End program
	return EXIT_SUCCESS;
}


/**
 This function processes the <priors> node of the XML problem definition file. The lower bound and upper
 bound information for each parameter are written in the lb and ub arrays. The arrays need to be preallocated
 to the correct dimension.
 */
int processPriorData (xmlDocPtr doc, double lb[], double ub[])
{
	// Initialize a new XPath evaluation context
	xmlXPathContextPtr ctxt = xmlXPathNewContext (doc);

	// Initialize a xmlXPathObjectPtr object
	xmlXPathObjectPtr result;

	// Apply XPath to find the set of <prior> nodes
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/priors/prior", ctxt);

	// Find out how many nodes were found (if any)
	size_t size = (result->nodesetval) ? result->nodesetval->nodeNr : 0;

	// Check that we have a non-empty nodeset
	if (size != R)
	{
		fprintf (stderr, "error in processPriorData: size mismatch\n");
		return GSL_FAILURE;
	}

	// Iterate through the resulting node set
	for (size_t i = 0; i < size; i++)
	{
		// Get first element child of current node -- this is the node <index>
		xmlNode * cur = xmlFirstElementChild (result->nodesetval->nodeTab[i]);

		// Extract the index of current parameter
		size_t j = (size_t) xmlXPathCastNodeToNumber (cur);

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// This is the node <dist>, which is currently not used.

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// Extract the lower bound encoded in the node <lb>
		double lbc = xmlXPathCastNodeToNumber (cur);

		// Advance to the next element sibling of the current node
		cur = xmlNextElementSibling (cur);

		// Extract the upper bound encoded in the node <ub>
		double ubc = xmlXPathCastNodeToNumber (cur);

		// Put stuff in place
		lb[j-1] = lbc;
		ub[j-1] = ubc;
	}

	// Free XML resources
	xmlXPathFreeObject (result);
	xmlXPathFreeContext (ctxt);

	return GSL_SUCCESS;
}


/**
 This function processes the <data> node of the XML problem definition file and reads
 the specified FCS files.
 */
int processDataSet (xmlDocPtr doc, gsl_vector * times, gsl_matrix * data)
{
	// Initialize a new XPath evaluation context
	xmlXPathContextPtr ctxt = xmlXPathNewContext (doc);

	// Initialize a xmlXPathObjectPtr object
	xmlXPathObjectPtr result;

	// Apply XPath to find the set of <timept> nodes
	result = xmlXPathEvalExpression ((xmlChar *) "/problem/data/timept", ctxt);

	// Find out how many nodes were found (if any)
	size_t size = (result->nodesetval) ? result->nodesetval->nodeNr : 0;

	// Check that we have a non-empty nodeset
	if (size != K)
	{
		fprintf (stderr, "error in processDataSet: size mismatch in time points\n");
		return GSL_FAILURE;
	}

	// Iterate through the resulting node set
	for (size_t i = 0; i < size; i++)
	{
		// Get first element child of current node -- this is the node <index>
		xmlNode * cur = xmlFirstElementChild (result->nodesetval->nodeTab[i]);
		// Extract the index of current time point
		size_t j = (size_t) xmlXPathCastNodeToNumber (cur);

		// Advance to the next element sibling of the current node -- this is <time>
		cur = xmlNextElementSibling (cur);
		// Extract the time of the measurement
		double t = xmlXPathCastNodeToNumber (cur);
		// Put the value in the appropriate place of the vector times
		gsl_vector_set (times, j-1, t);

		// Advance to the next element sibling of the current node -- this is <unit> -- currently unused
		cur = xmlNextElementSibling (cur);

		// Advance to the next element sibling of the current node -- this is <fcsfile>
		cur = xmlNextElementSibling (cur);
		// Extract the path to the FCS file containing the data
		xmlChar * path = xmlXPathCastNodeToString (cur);
		// Open the FCS file
		FILE * fcs = fopen (path, "rb");
		if (fcs == NULL)
		{
			fprintf (stderr, "error in processDataSet: could not open FCS file %d\n", (int) j);
			return GSL_FAILURE;
		}
		// Read the FCS file header
		float ver;
		long int dos;
		fscanf (fcs, "FCS%f    %*d%*d%ld%*d%*d%*d", &ver, &dos);

		// Initialize a new xmlXPathObjectPtr object
		xmlXPathObjectPtr result2;
		// Create the XPath expression to retrieve the node set of <measurement> for the current time point
		char xpath[60];
		sprintf (xpath, "/problem/data/timept[%d]/measuredspecies/measurement", (int) i+1);
		// Apply the XPath expression
		result2 = xmlXPathEvalExpression ((xmlChar *) xpath, ctxt);
		// Find out how many nodes were found (if any)
		size_t size2 = (result2->nodesetval) ? result2->nodesetval->nodeNr : 0;
		// Check that we have a non-empty nodeset
		if (size2 != P)
		{
			fprintf (stderr, "error in processDataSet: size mismatch in measured species\n");
			return GSL_FAILURE;
		}

		// Initialize arrays to contain indices and FCS parameter of the measured species
		size_t l[size2], p[size2];

		// Iterate through the resulting node set
		for (size_t k = 0; k < size2 ; k++)
		{
			// Get first element child of current node -- this is the node <index>
			xmlNode * cur2 = xmlFirstElementChild (result2->nodesetval->nodeTab[k]);
			// Extract the index of current measured species
			l[k] = (size_t) xmlXPathCastNodeToNumber (cur2);

			// Advance to the next element sibling of the current node -- this is <description> -- currently unused
			cur2 = xmlNextElementSibling (cur2);

			// Advance to the next element sibling of the current node -- this is <fcsparam>
			cur2 = xmlNextElementSibling (cur2);
			// Extract the parameter of current measured species
			p[k] = (size_t) xmlXPathCastNodeToNumber (cur2);
		}

		// Find out how many total events are stored in the FCS file
		int tot = 0;
		if (ver < 3.0)
			tot = fcs2_read_int_kw (fcs, "$TOT");
		else
			tot = fcs3_read_int_kw (fcs, "$TOT");
		// Check that there are enough events
		if (tot < M)
		{
			fprintf (stderr, "error in processDataSet: the FCS file %d does not contain enough events\n", (int) i+1);
			return GSL_FAILURE;
		}

		// Find out how many parameters are there in each event
		int par = 0;
		if (ver < 3.0)
			par = fcs2_read_int_kw (fcs, "$PAR");
		else
			par = fcs3_read_int_kw (fcs, "$PAR");

		// Position the file cursor at the beginning of the data segment
		fseek (fcs, dos, SEEK_SET);
		// Initialize read buffer
		float fbuf[par];
		// Read the first M events
		for (size_t m = 0; m < M; m++)
		{
			// Get a chunk of floats
			fread_floats_swap (fbuf, par, fcs);
			// Iterate through the parameters that have to be read
			for (size_t k = 0; k < size2; k++)
			{
				// Place the data in the proper spot
				gsl_matrix_set (data, m, P*(j-1)+l[k]-1, fbuf[p[k]]);
			}
		}

		// Close the FCS file
		fclose (fcs);

		// Free XML resources
		xmlXPathFreeObject (result2);
		xmlFree (path);
	}

	// Free XML resources
	xmlXPathFreeObject (result);
	xmlXPathFreeContext (ctxt);

	// Return
	return GSL_SUCCESS;
}
