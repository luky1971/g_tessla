/*
 * Copyright 2016 Ahnaf Siddiqui
 *
 * This program uses the GROMACS molecular simulation package API.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed at http://www.gromacs.org.
 */

#ifndef GKUT_LOG_H
#define GKUT_LOG_H

void init_log(const char *logfile, int argc, char *argv[]);
/* Opens the logfile and logs program command and initial time/date. Remember to close_log() at end of program.
 */

void close_log();
/* Closes the logfile.
 */

void print_log(const char *fmt, ...);
/* Prints to both stdout and the logfile.
 * Must call init_log at some point before calling this to print to logfile.
 */

void log_fatal(int fatal_errno, const char *file, int line, const char *fmt, ...);
/* Logs fatal error to logfile and also calls gmx_fatal
 * Hint: Use FARGS for first 3 arguments.
 * Must call init_log at some point before calling this to print to logfile.
 */

#endif // GKUT_LOG_H