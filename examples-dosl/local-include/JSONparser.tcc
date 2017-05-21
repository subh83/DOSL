/** **************************************************************************************
*                                                                                        *
*    A Ridiculously Simple JSON Parser for C++ (RSJP-cpp)                                *
*    Version 1.0a                                                                        *
*    ----------------------------------------------------------                          *
*    Copyright (C) 2016  Subhrajit Bhattacharya                                          *
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

#ifndef __DOSL_JSONPARSE_TCC
#define __DOSL_JSONPARSE_TCC

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>


char* JSONobjectbrackets = "{}";
char* JSONarraybrackets = "[]";
char JSONobjectassignment = ':';
char JSONarraydelimiter = ',';

std::vector<char*> JSONbrackets = {JSONobjectbrackets, JSONarraybrackets};
std::vector<char*> JSONcontainerquotes = {"\"\"", "''"};
char JSONcharescape = '\\';

// ============================================================

class JSONcontainer;
/* Use: JSONcontainer("JSON_string_data").as<JSONobject>()["keyName"].as<JSONarray>()[2].as<int>()
        JSONcontainer("JSON_string_data")["keyName"][2].as<int>()  */

// Helper preprocessor directives
#define jsonObject  as<JSONobject>()
#define jsonArray   as<JSONarray>()
#define jsonAs(t)   as<t>()

typedef std::map <std::string,JSONcontainer>    JSONobject;
typedef std::vector <JSONcontainer>             JSONarray;

class JSONcontainer {
public:
    std::string data; // can be object, vector or leaf data
    
    JSONcontainer () { }
    JSONcontainer (std::string str) : data (str) { }
    JSONcontainer (std::istream& is) {
        data = std::string ( (std::istreambuf_iterator<char>(is)), (std::istreambuf_iterator<char>()) );
    }
    
    template <class dataType>
    dataType as (void) { // specialized outside class declaration
        return dataType (data); // default behavior for unknown types: invoke 'dataType(std::string)'
    }
    
    // parsed data and opertor[]
    JSONobject parsed_obj;
    JSONcontainer& operator[] (std::string key);
    JSONarray parsed_array;
    JSONcontainer& operator[] (int indx);
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

int is_bracket (char c, std::vector<char*>& bracks, int indx=0) {
    for (int b=0; b<bracks.size(); ++b)
        if (c==bracks[b][indx]) 
            return (b);
    return (-1);
}

std::vector<std::string> split_JSON_array (std::string str) {
    // splits, while respecting brackets and escapes
    std::vector<std::string> ret;
    
    std::string current;
    std::vector<int> bracket_stack;
    std::vector<int> quote_stack;
    bool escape_active = false;
    int bi;
    
    for (int a=0; a<str.length(); ++a) { // *
        
        // delimiter
        if ( bracket_stack.size()==0  &&  quote_stack.size()==0  &&  str[a]==JSONarraydelimiter ) {
            ret.push_back (current);
            current.clear(); bracket_stack.clear(); quote_stack.clear(); escape_active = false;
            continue; // to *
        }
        
        // ------------------------------------
        // checks for string
        
        if (quote_stack.size() > 0) { // already inside string
            if (str[a]==JSONcharescape)  // an escape character
                escape_active = !escape_active;
            else if (!escape_active  &&  str[a]==JSONcontainerquotes[quote_stack.back()][1] ) { // close quote
                quote_stack.pop_back();
                escape_active = false;
            }
            else
                escape_active = false;
            
            current.push_back (str[a]);
            continue; // to *
        }
        
        if (quote_stack.size()==0) { // check for start of string
            if ((bi = is_bracket (str[a], JSONcontainerquotes)) >= 0) {
                quote_stack.push_back (bi);
                current.push_back (str[a]);
                continue; // to *
            }
        }
        
        // ------------------------------------
        // checks for brackets
        
        if ( bracket_stack.size()>0  &&  str[a]==JSONbrackets[bracket_stack.back()][1] ) { // check for closing bracket
            bracket_stack.pop_back();
            current.push_back (str[a]);
            continue;
        }
        
        if ((bi = is_bracket (str[a], JSONbrackets)) >= 0) {
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

// JSONobject
template <>
std::map <std::string,JSONcontainer>  JSONcontainer::as<JSONobject> (void) {
    std::map <std::string,JSONcontainer> ret;
    
    std::string content = strtrim (strtrim (data, " {", "l" ), " }", "r" );
    
    std::vector<std::string> nvPairs = split_JSON_array (content);
    for (int a=0; a<nvPairs.size(); ++a) {
        std::size_t assignmentPos = nvPairs[a].find (JSONobjectassignment);
        ret.insert (make_pair( 
                            strip_outer_quotes (nvPairs[a].substr (0,assignmentPos) ) ,
                            JSONcontainer (strtrim (nvPairs[a].substr (assignmentPos+1) ) )
                   ) );
    }
    
    return (ret);
}

JSONcontainer& JSONcontainer::operator[] (std::string key) {
    if (parsed_obj.empty())
        parsed_obj = as<JSONobject>();
    return (parsed_obj[key]);
}

// ------------------------------------

// JSONarray
template <>
std::vector <JSONcontainer>  JSONcontainer::as<JSONarray> (void) {
    std::vector <JSONcontainer> ret;
    
    std::string content = strtrim (strtrim (data, " [", "l" ), " ]", "r" );
    
    std::vector<std::string> nvPairs = split_JSON_array (content);
    for (int a=0; a<nvPairs.size(); ++a) 
        ret.push_back (JSONcontainer (strtrim (nvPairs[a]) ) );
    
    return (ret);
}

JSONcontainer& JSONcontainer::operator[] (int indx) {
    if (parsed_array.empty())
        parsed_array = as<JSONarray>();
    return (parsed_array[indx]);
}

// ------------------------------------

// String
template <>
std::string  JSONcontainer::as<std::string> (void) {
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
int  JSONcontainer::as<int> (void) {
    return (atoi (strip_outer_quotes(data).c_str() ) );
}

// double
template <>
double  JSONcontainer::as<double> (void) {
    return (atof (strip_outer_quotes(data).c_str() ) );
}

// bool
template <>
bool  JSONcontainer::as<bool> (void) {
    std::string cleanData = strip_outer_quotes (data);
    if (cleanData=="true" || cleanData=="TRUE" || cleanData=="True" || atoi(cleanData.c_str())!=0) return (true);
    return (false);
}

#endif
