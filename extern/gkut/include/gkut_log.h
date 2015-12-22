/*
 * Copyright Ahnaf Siddiqui
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