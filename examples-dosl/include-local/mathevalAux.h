/** **************************************************************************************
*                                                                                        *
*    Part of                                                                             *
*    Discrete Optimal search Library (DOSL)                                              *
*    A template-based C++ library for discrete search                                    *
*    Version 3.x                                                                         *
*    ----------------------------------------------------------                          *
*    Copyright (C) 2017  Subhrajit Bhattacharya                                          *
*                                                                                        *
*    This program is free software: you can redistribute it and/or modify                *
*    it under the terms of the GNU General Public License as published by                *
*    the Free Software Foundation, either version 3 of the License, or                   *
*    (at your option) any later version.                                                 *
*                                                                                        *
*    This program is distributed in the hope that it will be useful,                     *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of                      *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                       *
*    GNU General Public License for more details <http://www.gnu.org/licenses/>.         *
*                                                                                        *
*                                                                                        *
*    Contact:  subhrajit@gmail.com                                                       *
*              https://www.lehigh.edu/~sub216/ , http://subhrajit.net/                   *
*                                                                                        *
*                                                                                        *
*************************************************************************************** **/

#ifndef __DOSL_MATHEVAL_AUX_TCC
#define __DOSL_MATHEVAL_AUX_TCC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <matheval.h>
#include <csignal>


class MathEvaluator {
private:
    char* expression;
    char**  symbols;
    double* values;
    void* evaluator_obj;
    
    bool initiated;
    int n_vars;
    std::unordered_map<std::string,double*> varps;
    
    // -------------------------------
    // destructor
    
    void destroy() {
        if (initiated) {
            evaluator_destroy (evaluator_obj);
            delete[] expression;
            for (int a=0; a<n_vars; ++a)
                delete[] symbols[a];
            delete[] symbols;
            delete[] values;
            initiated = false;
        }
    }
    
public:
    
    MathEvaluator () : initiated (false) { }
    
    MathEvaluator (std::string expr, std::vector<std::string> varlist ) : initiated (false) {
        init (expr, varlist);
    }
    
    MathEvaluator (std::string expr, std::unordered_map<std::string,double> var_vals=std::unordered_map<std::string,double>() ) 
                : initiated (false) {
        init (expr, var_vals);
    }
    
    // -----------
    
    bool empty(void) { return (!initiated); }
    
    // -----------
    
    void init (std::string expr, std::unordered_map<std::string,double> var_vals=std::unordered_map<std::string,double>()) {
        std::vector<std::string> varlist;
        for (auto it=var_vals.begin(); it!=var_vals.end(); ++it)
            varlist.push_back(it->first);
        init (expr, varlist); // will call destroy()
        for (auto it=var_vals.begin(); it!=var_vals.end(); ++it)
            set_var (it->first, it->second);
    }
    
    void init (std::string expr, std::vector<std::string> varlist ) {
        destroy(); // if already initiated
        
        expression = new char[expr.size()+1];
        strcpy (expression, expr.c_str());
        evaluator_obj = evaluator_create (expression);
        if (!evaluator_obj) {
            printf("Error in creating evaluator_obj. Possible syntax error in expression: '%s'.\n", expr.c_str());
            std::abort();
        }
        
        n_vars = varlist.size();
        symbols = new char*[n_vars];
        values = new double[n_vars];
        for (int a=0; a<n_vars; ++a) {
            symbols[a] = new char[varlist[a].size()+1];
            strcpy (symbols[a], varlist[a].c_str());
            varps[varlist[a]] = (values+a);
        }
        
        initiated = true;
    }
    
    // -----------
    
    void set_var (std::string sym, double val) {
        *(varps[sym]) = val;
    }
    
    double eval (std::unordered_map<std::string,double> var_vals=std::unordered_map<std::string,double>() ) {
        for (auto it=var_vals.begin(); it!=var_vals.end(); ++it)
            set_var (it->first, it->second);
        // debug:
        /*printf("Evaluating '%s' with... ", evaluator_get_string(evaluator_obj));
        for (int a=0; a<n_vars; ++a)
            printf ("%s = %f; ", symbols[a], values[a]);
        printf("\n");*/
        
        return (evaluator_evaluate (evaluator_obj, n_vars, symbols, values));
    }
    
    // -------------------------------
    // copy: TODO
    
    MathEvaluator (const MathEvaluator& other) {
        if (other.initiated) {
            std::unordered_map<std::string,double> var_vals;
            for (auto it=other.varps.begin(); it!=other.varps.end(); ++it)
                var_vals[it->first] = *(it->second);
            init (other.expression, var_vals);
        }
        else
            initiated = false;
    }
    
    MathEvaluator& operator= (const MathEvaluator& other) {
        if (other.initiated) {
            std::unordered_map<std::string,double> var_vals;
            for (auto it=other.varps.begin(); it!=other.varps.end(); ++it)
                var_vals[it->first] = *(it->second);
            init (other.expression, var_vals);
        }
        else
            initiated = false;
        return *this;
    }
    
    // -------------------------------
    // destructor
    
    ~MathEvaluator() {
        destroy();
    }

};

#endif
