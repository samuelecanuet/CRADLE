#ifndef MESSENGER
#define MESSENGER

#include <string>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <dirent.h>
#include "TLatex.h"

// TERMINAL COLORS //
inline std::string RED = "\033[1;31m";
inline std::string GREEN = "\033[1;32m";
inline std::string YELLOW = "\033[1;33m";
inline std::string BLUE = "\033[1;34m";
inline std::string MAGENTA = "\033[1;35m";
inline std::string CYAN = "\033[1;36m";
inline std::string WHITE = "\033[1;37m";
inline std::string RESET = "\033[0m";

inline int GENERAL_INDENT = 12;

// GENERAL MESSENGER FUNCTIONS //

inline void Error(const std::string &message)
{
    std::cout << RED << std::left << std::setw(GENERAL_INDENT) << " <ERROR>" << message << RESET << std::endl;
    exit(0);
}

inline void Error(const char *message)
{
    Error(std::string(message));
}

inline void Warning(const std::string &message, int indent = 0)
{
    std::cout << YELLOW << std::left << std::setw(GENERAL_INDENT) << " <WARNING>";
    if (indent > 0)
    {
        for (int i = 0; i < indent; i++)
            std::cout << " ";
        std::cout << "├─";
    }
    std::cout << message << RESET << std::endl;
}

inline void Warning(const char *message, int indent = 0)
{
    Warning(std::string(message), indent);
}

inline void Info(const std::string &message, int indent = 0)
{
    std::cout << BLUE << std::left << std::setw(GENERAL_INDENT) << " <INFO>  ";
    if (indent > 0)
    {
        for (int i = 0; i < indent; i++)
            std::cout << " ";
        std::cout << "├─";
    }
    std::cout << message << RESET << std::endl;
}

inline void Info(const char *message, int indent = 0)
{
    Info(std::string(message), indent);
}

inline void Success(const std::string &message, int indent = 0)
{
    std::cout << GREEN << std::left << std::setw(GENERAL_INDENT) << " <SUCCESS>";
    if (indent > 0)
    {
        for (int i = 0; i < indent; i++)
            std::cout << " ";
        std::cout << "├─";
    }
    std::cout << message << RESET << std::endl;
}

inline void Success(const char *message, int indent = 0)
{
    Success(std::string(message), indent);
}

inline void Start(const std::string &message, int indent = 0)
{
    std::cout << CYAN << std::left << std::setw(GENERAL_INDENT) << " <START>  ";
    if (indent > 0)
    {
        for (int i = 0; i < indent; i++)
            std::cout << " ";
        std::cout << "├─";
    }
    std::cout << message << RESET << std::endl;
}

inline void Start(const char *message, int indent = 0)
{
    Start(std::string(message), indent);
}

// PERSONNALIZED MESSENGER FUNCTIONS //
inline std::string GetColorWithstring(std::string color)
{
    // converting to lower case
    std::transform(color.begin(), color.end(), color.begin(), [](unsigned char c) { return std::tolower(c); });
    if (color == "red")
        return RED;
    else if (color == "green")
        return GREEN;
    else if (color == "yellow")
        return YELLOW;
    else if (color == "blue")
        return BLUE;
    else if (color == "magenta")
        return MAGENTA;
    else if (color == "cyan")
        return CYAN;
    else if (color == "white")
        return WHITE;
    else
        return RESET;
}

inline void Message(std::string type, const std::string &message, int indent = 0, std::string color = "")
{
    color = GetColorWithstring(color);
    std::cout << color;
    if (type != "")
    {
        std::ostringstream typeSS;
        typeSS << " <" << type << "> ";
        std::cout << std::left << std::setw(GENERAL_INDENT) << typeSS.str();
    }

    if (indent > 0)
    {
        for (int i = 0; i < indent; i++)
            std::cout << " ";
        std::cout << "├─";
    }
    std::cout << message << RESET << std::endl;
}

inline void ProgressBar(int cEntry, int TotalEntries, clock_t start, std::string Prefix = "", int Step=10000, int threads=1)
{
  if (cEntry % Step == 0)
  {
    clock_t Current = clock();
    double Frac = 1.0 * cEntry / TotalEntries;
    double Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC / (double)threads);
    double TimeLeft = Timeclock * (1. / Frac - 1.);

    std::cout << ("\r \u23F3 "+Prefix).c_str()
         << " - "
         << Form("%4.1f", 100. * cEntry / TotalEntries) << " %"
         << " - "
         << " Time Left : " << Form("%2d min ", (int)TimeLeft / 60)
         << Form("%02d sec", (int)TimeLeft % 60)
         << std::flush;
  }

  if (cEntry == TotalEntries)
  {
    std::cout << std::endl;
  }  
}

#endif // MESSENGER