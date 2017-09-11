/** **************************************************************************************
*                                                                                        *
*    A Ridiculously Simple JSON parser for C++ (RSJp-cpp)                                *
*    Version 1.0b                                                                        *
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

#ifndef __DOSL_RSJPARSE_TCC
#define __DOSL_RSJPARSE_TCC

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>


char const* RSJobjectbrackets = "{}";
char const* RSJarraybrackets = "[]";
char RSJobjectassignment = ':';
char RSJarraydelimiter = ',';

std::vector<char const*> RSJbrackets = {RSJobjectbrackets, RSJarraybrackets};
std::vector<char const*> RSJstringquotes = {"\"\"", "''"};
char RSJcharescape = '\\';

std::string RSJlinecommentstart = "//";

// ============================================================

class RSJresource;
/* Use: RSJresource("RSJ_string_data").as<RSJobject>()["keyName"].as<RSJarray>()[2].as<int>()
        RSJresource("RSJ_string_data")["keyName"][2].as<int>()  */

// Helper preprocessor directives
#define rsjObject  as<RSJobject>()
#define rsjArray   as<RSJarray>()
#define rsjAs(t)   as<t>()

typedef std::map <std::string,RSJresource>    RSJobject;
typedef std::vector <RSJresource>             RSJarray;

class RSJresource {
private:
    // main data
    std::string data; // can be object, vector or leaf data
    bool _exists;      // whether the RSJ resource exists.
    
    // parsed data
    RSJobject parsed_obj;
    RSJarray parsed_array;
    
public:
    // constructor
    RSJresource () : _exists (false) { } // no data field.
    
    RSJresource (std::string str) : data (str), _exists (true) {
        if (data.empty()) _exists = false;
        else _exists = true;
    }
    
    RSJresource (std::istream& is) {
        data = std::string ( (std::istreambuf_iterator<char>(is)), (std::istreambuf_iterator<char>()) );
        if (data.empty()) _exists = false;
        else _exists = true;
    }
    
    // access raw data
    std::string& raw_data (void) { return (data); }
    bool exists (void) { return (_exists); }
    
    // as
    template <class dataType>
    dataType as (const dataType& def = dataType()) { // specialized outside class declaration
        if (!exists()) return (def);
        return dataType (data); // default behavior for unknown types: invoke 'dataType(std::string)'
    }
    
    // opertor[]
    RSJresource& operator[] (std::string key); // object
    RSJresource& operator[] (int indx); // array
};

// ============================================================
// Direct string manipulation functions

std::string strtrim (std::string str, std::string chars=" \n\r", std::string opts="lr") {
    if (str.empty()) return(str);
    
    if (opts.find('l')!=std::string::npos) { // left trim
        int p;
        for (p=0; p<str.length(); ++p)
            if (chars.find(str[p])==std::string::npos) break;
        str.erase (0, p);
    }
    
    if (opts.find('r')!=std::string::npos) { // right trim
        int q, strlenm1=str.length()-1;
        for (q=0; q<str.length(); ++q)
            if (chars.find(str[strlenm1-q])==std::string::npos) break;
        str.erase (str.length()-q, q);
    }
    
    return (str);
}

std::string strip_outer_quotes (std::string str, char* qq=NULL) {
    str = strtrim (str);
    
    std::string ret = strtrim (str, "\"");
    if (ret==str) {
        ret = strtrim (str, "'");
        if (qq && ret!=str) *qq = '\'';
    }
    else if (qq)
        *qq = '"';
    
    return (ret);
}

// ----------------

int is_bracket (char c, std::vector<char const*>& bracks, int indx=0) {
    for (int b=0; b<bracks.size(); ++b)
        if (c==bracks[b][indx]) 
            return (b);
    return (-1);
}

std::vector<std::string> split_RSJ_array (std::string str) {
    // splits, while respecting brackets and escapes
    std::vector<std::string> ret;
    
    std::string current;
    std::vector<int> bracket_stack;
    std::vector<int> quote_stack;
    bool escape_active = false;
    int bi;
    
    for (int a=0; a<str.length(); ++a) { // *
        
        // delimiter
        if ( bracket_stack.size()==0  &&  quote_stack.size()==0  &&  str[a]==RSJarraydelimiter ) {
            ret.push_back (current);
            current.clear(); bracket_stack.clear(); quote_stack.clear(); escape_active = false;
            continue; // to *
        }
        
        // ------------------------------------
        // checks for string
        
        if (quote_stack.size() > 0) { // already inside string
            if (str[a]==RSJcharescape)  // an escape character
                escape_active = !escape_active;
            else if (!escape_active  &&  str[a]==RSJstringquotes[quote_stack.back()][1] ) { // close quote
                quote_stack.pop_back();
                escape_active = false;
            }
            else
                escape_active = false;
            
            current.push_back (str[a]);
            continue; // to *
        }
        
        if (quote_stack.size()==0) { // check for start of string
            if ((bi = is_bracket (str[a], RSJstringquotes)) >= 0) {
                quote_stack.push_back (bi);
                current.push_back (str[a]);
                continue; // to *
            }
        }
        
        // ------------------------------------
        // checks for comments
        
        if (quote_stack.size()==0) { // comment cannot start inside string
            
            // single-line commenst
            if (str.compare (a, RSJlinecommentstart.length(), RSJlinecommentstart) == 0) {
                // ignore until end of line
                int newline_pos = str.find ("\n", a);
                if (newline_pos == std::string::npos)
                    newline_pos = str.find ("\r", a);
                
                if (newline_pos != std::string::npos)
                    a = newline_pos; // point to the newline character (a will be incremented)
                else // the comment continues until EOF
                    a = str.length();
                continue;
            }
        }
        
        // ------------------------------------
        // checks for brackets
        
        if ( bracket_stack.size()>0  &&  str[a]==RSJbrackets[bracket_stack.back()][1] ) { // check for closing bracket
            bracket_stack.pop_back();
            current.push_back (str[a]);
            continue;
        }
        
        if ((bi = is_bracket (str[a], RSJbrackets)) >= 0) {
            bracket_stack.push_back (bi);
            current.push_back (str[a]);
            continue; // to *
        }
        
        // ------------------------------------
        // otherwise
        current.push_back (str[a]);
    }
    
    if (current.length() > 0)
        ret.push_back (current);
    
    return (ret);
}


// ============================================================
// Specialized .as() member functions

// RSJobject
template <>
std::map <std::string,RSJresource>  RSJresource::as<RSJobject> (const RSJobject& def) {
    if (!exists()) return (def);
    
    std::map <std::string,RSJresource> ret;
    std::string content = strtrim (strtrim (strtrim(data), " {", "l" ), " }", "r" );
    
    std::vector<std::string> nvPairs = split_RSJ_array (content);
    for (int a=0; a<nvPairs.size(); ++a) {
        std::size_t assignmentPos = nvPairs[a].find (RSJobjectassignment);
        ret.insert (make_pair( 
                            strip_outer_quotes (nvPairs[a].substr (0,assignmentPos) ) ,
                            RSJresource (strtrim (nvPairs[a].substr (assignmentPos+1) ) )
                   ) );
    }
    
    return (ret);
}

RSJresource& RSJresource::operator[] (std::string key) {
    if (parsed_obj.empty())
        parsed_obj = as<RSJobject>();
    return (parsed_obj[key]); // will return empty resource if key does not exist
}

// ------------------------------------

// RSJarray
template <>
std::vector <RSJresource>  RSJresource::as<RSJarray> (const RSJarray& def) {
    if (!exists()) return (def);
    
    std::vector <RSJresource> ret;
    std::string content = strtrim (strtrim (strtrim(data), " [", "l" ), " ]", "r" );
    
    std::vector<std::string> nvPairs = split_RSJ_array (content);
    for (int a=0; a<nvPairs.size(); ++a) 
        ret.push_back (RSJresource (strtrim (nvPairs[a]) ) );
    
    return (ret);
}

RSJresource& RSJresource::operator[] (int indx) {
    if (parsed_array.empty())
        parsed_array = as<RSJarray>();
    if (indx>=parsed_array.size())
        parsed_array.resize(indx+1); // insert empty resources
    return (parsed_array[indx]);
}

// ------------------------------------

// String
template <>
std::string  RSJresource::as<std::string> (const std::string& def) {
    if (!exists()) return (def);
    
    char qq = '\0';
    std::string ret = strip_outer_quotes (data, &qq);
    
    std::vector< std::vector<std::string> > escapes = { {"\\n","\n"}, {"\\r","\r"}, {"\\t","\t"}, {"\\\\","\\"} };
    if (qq=='"')
        escapes.push_back ({"\\\"","\""});
    else if (qq=='\'')
        escapes.push_back ({"\\'","'"});
    
    for (int a=0; a<escapes.size(); ++a)
        for ( std::size_t start_pos=ret.find(escapes[a][0]); start_pos!=std::string::npos; start_pos=ret.find(escapes[a][0],start_pos) ) {
            ret.replace (start_pos, escapes[a][0].length(), escapes[a][1]);
            start_pos += escapes[a][1].length();
        }
    
    return (ret);
}

// integer
template <>
int  RSJresource::as<int> (const int& def) {
    if (!exists()) return (def);
    return (atoi (strip_outer_quotes(data).c_str() ) );
}

// double
template <>
double  RSJresource::as<double> (const double& def) {
    if (!exists()) return (def);
    return (atof (strip_outer_quotes(data).c_str() ) );
}

// bool
template <>
bool  RSJresource::as<bool> (const bool& def) {
    if (!exists()) return (def);
    std::string cleanData = strip_outer_quotes (data);
    if (cleanData=="true" || cleanData=="TRUE" || cleanData=="True" || atoi(cleanData.c_str())!=0) return (true);
    return (false);
}

#endif
