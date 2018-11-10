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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

// SGL header:
#include <sgl/sgl>

int main(int argc, char *argv[])
{
    // ===============================================================================
    // Example 0:
    
    std::cout << "\n---------\nEXAMPLE 0: " << std::endl;
    std::cout << RSJresource("{'some_name': string_data, keyName: [2,3,5,7]}")["keyName"][2].as<int>() << std::endl;
    
    
    // ===============================================================================
    // Example 1:
    
    std::string    str = "{'animal':cat, coordinates: [2, 5, 8], height: 1, \nis_vicious: false, comment:'It\\'s in fact quite...\\t adorable.' }";
    RSJresource    my_resource (str); // RSJ parser.
    
    std::cout << "\n---------\nEXAMPLE 1: " << std::endl;
    std::cout << "The JSON string:\n\n" << str << "\n" << std::endl;
    
    std::cout << "Iterating over the object fields:" << std::endl;
    RSJobject obj_map = my_resource.as<RSJobject>(); // RSJobject = std::map <std::string,RSJresource>
    for (auto it=obj_map.begin(); it!=obj_map.end(); ++it)
        std::cout << "\t" << it->first << " => " << it->second.raw_data() << std::endl;
    std::cout << std::endl;
    
    std::cout << "Some specific queries:" << std::endl;
    std::cout << "\tThe animal is: " << my_resource["animal"].as<std::string>("dog") << std::endl;
    std::cout << "\tIts Y coordinate is: " << my_resource["coordinates"][1].as<int>() << std::endl;
    std::cout << "\tIts Z coordinate is: " << my_resource["coordinates"][2].as<double>() << std::endl;
    std::cout << "\tIs it vicious? " << my_resource["is_vicious"].as<bool>() << std::endl;
    std::cout << "\tComment: " << my_resource["comment"].as<std::string>() << std::endl;
    
    
    if (my_resource["length"].exists())
        std::cout << "\tLength: " << my_resource["length"].as<int>() << std::endl;
    else 
        std::cout << "\tLength: [does not exist.]"  << std::endl;
    
    if (my_resource["height"].exists())
        std::cout << "\tHeight: " << my_resource["height"].as<int>() << std::endl;
    else
        std::cout << "\tHeight: [does not exist.]"  << std::endl;
    
    int default_width = -1;
    std::cout << "\tWidth: " << my_resource["width"].as<int>(default_width) << std::endl;
}


