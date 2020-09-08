#pragma once
#include "hdf_wrapper.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>


auto readDataSetREAL4=[](auto &gr,auto &array,std::string name,size_t N)
{
    std::cout<<"Reading from file REAL4 size="<<N<<" array name="<<name<<"\n";
    for(size_t ii=0;ii<gr.size();ii++)
        if(gr.get_link_name(ii).compare(name)==0){
            auto dt = gr.open_dataset(name.c_str());
            std::vector<double> tempArray(N*4);
            dt.read(&tempArray[0]);
            for(size_t i=0;i<N;i++)
            {
                array[i]=REAL4({tempArray[i*4],tempArray[i*4+1],tempArray[i*4+2],tempArray[i*4+3]});
            }
        }
};
auto readDataSetREAL=[](auto &gr,auto &array,std::string name,size_t N)
{
    std::cout<<"Reading from file REAL size="<<N<<" array name="<<name<<"\n";
    for(size_t ii=0;ii<gr.size();ii++)
        if(gr.get_link_name(ii).compare(name)==0){
            auto dt = gr.open_dataset(name.c_str());
            std::vector<double> tempArray(N);
            dt.read(&tempArray[0]);
            for(size_t i=0;i<N;i++)
            {
                array[i]=REAL(tempArray[i]);
            }
        }
};
auto readDataSetINT2=[](auto &gr,auto &array,std::string name,size_t N)
{
    std::cout<<"Reading from file INT2 size="<<N<<" array name="<<name<<"\n";
    for(size_t ii=0;ii<gr.size();ii++)
        if(gr.get_link_name(ii).compare(name)==0){
            auto dt = gr.open_dataset(name.c_str());
            std::vector<int> tempArray(N*2);
            dt.read(&tempArray[0]);
            for(size_t i=0;i<N;i++)
            {
                array[i]=INT2({tempArray[i*2],tempArray[i*2+1]});
            }
        }
};
auto readDataSetINT4=[](auto &gr,auto &array,std::string name,size_t N)
{
    std::cout<<"Reading from file INT4 size="<<N<<" array name="<<name<<"\n";
    for(size_t ii=0;ii<gr.size();ii++)
        if(gr.get_link_name(ii).compare(name)==0){
            auto dt = gr.open_dataset(name.c_str());
            std::vector<int> tempArray(N*4);
            dt.read(&tempArray[0]);
            for(size_t i=0;i<N;i++)
            {
                array[i]=INT4({tempArray[i*4],tempArray[i*4+1],tempArray[i*4+2],tempArray[i*4+3]});
            }
        }
};
auto readDataSetINT=[](auto &gr,auto &array,std::string name,size_t N)
{
    std::cout<<"Reading from file INT size="<<N<<" array name="<<name<<"\n";
    for(size_t ii=0;ii<gr.size();ii++)
        if(gr.get_link_name(ii).compare(name)==0){
            auto dt = gr.open_dataset(name.c_str());
            std::vector<int> tempArray(N);
            dt.read(&tempArray[0]);
            for(size_t i=0;i<N;i++)
            {
                array[i]=tempArray[i];
            }
        }
};


auto readDataSetUCHAR4=[](auto &gr,auto &array,std::string name,size_t N)
{
    std::cout<<"Reading from file UCHAR4 size="<<N<<" array name="<<name<<"\n";
    for(size_t ii=0;ii<gr.size();ii++)
        if(gr.get_link_name(ii).compare(name)==0){
            auto dt = gr.open_dataset(name.c_str());
            std::vector<UCHAR> tempArray(N*4);
            dt.read(&tempArray[0]);
            for(size_t i=0;i<N;i++)
            {
                array[i]=UCHAR4({tempArray[i*4],tempArray[i*4+1],tempArray[i*4+2],tempArray[i*4+3]});
            }
        }
};
auto readDataSetUCHAR=[](auto &gr,auto &array,std::string name,size_t N)
{
    std::cout<<"Reading from file UCHAR size="<<N<<" array name="<<name<<"\n";
    for(size_t ii=0;ii<gr.size();ii++)
        if(gr.get_link_name(ii).compare(name)==0){
            auto dt = gr.open_dataset(name.c_str());
            std::vector<UCHAR> tempArray(N);
            dt.read(&tempArray[0]);
            for(size_t i=0;i<N;i++)
            {
                array[i]=tempArray[i];
            }
        }
};


