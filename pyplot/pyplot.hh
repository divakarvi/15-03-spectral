#ifndef __PYPLOT19JULY2013__
#define __PYPLOT19JULY2013__
#include <cstdio>
#include <unistd.h>

enum pipe_type {PLTOFF, PLTON};

const int MAX_NUM_PYPLT_CMDS = 400;
const int MAX_CMD_PYPLT_LEN = 3000;

class PyPlot{
private:
	FILE *pypipe; //pipe to Python
	/*
	 * state = 0 (no line drawn)
         *       = 1 (line drawn, but figure not saved/shown)
         *       = 2 (figure saved/shown, should possibly save script and exit)
	 */
	int state;
	/*
	 * if nonzero no pipe is opened to python
	 */
	enum pipe_type pipe_state;
	/*
	 * figure will be saved in FIGS/"name.pdf"
	 * if asked to do so, data & script will be saved in FIGS/name_*
	 * normally data is saved in FIGS/name_* and deleted
	 */
	char name[30];
	/*
	 * array of commands issued to python
	 * at most MAX_NUM commands are allowd
	 */
	char *cmd[MAX_NUM_PYPLT_CMDS];
	char cmdstr[MAX_CMD_PYPLT_LEN];
	/*
	 * cmdnum = 0 initially
	 *        = index of current/next command
	 * incremented after cmd is issued
	 */
	int cmdnum;
	/*
	 * linenum = index of current line
	 * appended to name_ to save data for Python
	 * all saved data is removed unless specified otherwise
	 */
	int linenum;
	/*
	 * savedata = 0 then line data is removed by destructor
	 *          = 1 then line data is not removed from FIGS/name_*
	 */
	int savedata;
	/*
	 * called by plot(), semilogx/y, loglog
	 */
	void prep_line(double *x, double *y, int n);
	void issue_command(const char *cstr);
public:
	/*
	 * name must be less than 25 chars
	 * pipe == PLTOFF then mpl.backend is set to PDF
	 */
	PyPlot(const char *namei, enum pipe_type pipe=PLTON);
	~PyPlot();
	/*
	 * functions for drawing lines and setting their properties
	 */
	void plot(double *x, double *y, int n);
	void plot(double *y, int n);
	void semilogx(double *x, double *y, int n);
	void semilogy(double *x, double *y, int n);
	void semilogy(double *y, int n);
	void loglog(double *x, double *y, int n);
	void linestyle(const char* s);
	void linewidth(const char* s);
	void linecolor(const char* s);
	void markertype(const char* s); 
	void markercolor(const char* s);
	void markersize(const char* s);
	/*
	 * functions for specifying axes and their properties
	 */
	void axis(); //"tight"
	void axis(double x0, double x1, double y0, double y1);
	void title(const char* s);
	void xticks(double *ticks, int n);
	void yticks(double *ticks, int n);
	void ticksize(const char *s);
	/*
	 * issues command to python pipe verbatim
	 * last char of s must be '\n'
	 * for python syntax inside s, ax and l are always available as
	 * current axis and current line
	 */
	void pycmd(const char *s);
	/*
	 * function for showing/output
	 * if PLTON then show plot
	 * if PLTOFF then output to pdf
	 */
	void show();

	/*
	 * save python script in FIGS/
	 * must be called AFTER show() at the very end
	 */
	void savescript();
};
#endif
