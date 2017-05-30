#include <time.h>
#include <ctime>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include "vtkNew.h"


#include "mpi.h"










#include <sstream>
#include <iostream>
#include <sstream>
#include <ctype.h>
#include <string.h>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>


#include <ctype.h>
#include <string.h>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <utility>


#include <sys/stat.h>
#include <vector>      // for NodeDataInfoTemp & CellDataInfoTemp
#include <map>         // for TimeStepValuesMap & NumberOfPolygonsMap
#include <string>      // for TimeStepValuesMap & NumberOfPolygonsMap
 
 
 
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator> 

#include <cstring>
#include <string>







/*
extern"C" {
   void dodo2vtk_Init(char*, int, int*, int*, int*);
   void dodo2vtk_GetElementVolume(int num, double *volume);
   void dodo2vtk_GetElementType(int num, int *vtktype);
   void dodo2vtk_GetElementKonnec(int num, int vtkkonec[]);
   void dodo2vtk_GetNodeValueVector(int, char*, int, double vect[3]);
   void dodo2vtk_GetElementValueTensor2s(int num,char* cval,int lcval,double tens[6]);
   void dodo2vtk_GetElementMaterialValueScalar(int num,char* cval,int lcval,char* cmat,int lcmat,double* scal);
   void dodo2vtk_GetElementMaterialValueTensor1(int num,char* cval,int lcval,char* cmat,int lcmat,double tens[3]);
   void dodo2vtk_ReadNextFrame(int* step, int* cyc, int* inc, double* time, double* time_step, double* time_cycle, int* ierr);
   void dodo2vtk_End();
   void dodo2vtk_Rewind();
}
*/

#include "dodoF2C.hpp"
#include "cal.h"





int GetMesh(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nnodes, int nelts)
{

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    std::string name="COO";
    char *cname = new char[name.length() + 1];
    std::strcpy(cname, name.c_str());

    //std::cout << name.length() << " " << cname << std::endl;

    for (int i=0; i<nnodes; i++)
    {
     double vect[3]={0,0,0};
     dodoF2C_GetNodeValueVector(i+1,cname,name.length(),vect);
     //std::cout << vx << " " << vy << " " << vz << std::endl;
     points->InsertNextPoint(vect);
    }
    unstructuredGrid->SetPoints(points);
    delete[] cname;

    int etype;
    int konnec[100]; 
    for (int i=0; i<nelts;i++)
    {
     dodoF2C_GetElementTypeVTK(i+1,&etype);
     dodoF2C_GetElementKonnecVTK(i+1,konnec);

     switch (etype)
     {
      case 5: {
       vtkIdType ptIds[] = {konnec[0]-1, konnec[1]-1, konnec[2]-1};
       unstructuredGrid->InsertNextCell( VTK_TRIANGLE, 3, ptIds );
       break;
      }
      case 10: {
       vtkIdType ptIds[] = {konnec[0]-1, konnec[1]-1, konnec[2]-1, konnec[3]-1};
       unstructuredGrid->InsertNextCell( VTK_TETRA, 4, ptIds );
       break;
      }
      case 12: {
       vtkIdType ptIds[] = {konnec[0]-1, konnec[1]-1, konnec[2]-1, konnec[3]-1, konnec[4]-1, konnec[5]-1, konnec[6]-1, konnec[7]-1};
       unstructuredGrid->InsertNextCell( VTK_HEXAHEDRON, 8, ptIds );
       break;
      }
      default: {
       std::cerr << "ERROR, ELEMENT TYPE " << etype << "NOT YET IMPLEMENTED";
      }
     }
    }



  return 0;
}


int GetNodeValueVector(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nnodes, std::string name, std::string newname)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());
    {
     vtkSmartPointer<vtkDoubleArray> vectors =vtkSmartPointer<vtkDoubleArray>::New();
     vectors->SetNumberOfComponents(3);
     vectors->SetName(newname.c_str());
     for (int i=0;i<nnodes;i++)
     {
      double vect[3]={0,0,0};
      dodoF2C_GetNodeValueVector(i+1,cname,name.size(),vect);
      vectors->InsertNextTupleValue(vect);
     }
     unstructuredGrid->GetPointData()->AddArray(vectors);
    }
    delete[] cname;

  return 0;
}


int GetElementValueTensor2s(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts, std::string name, std::string newname)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());

    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(6); //we will have only 1 value associated with the triangle
    vals->SetName(newname.c_str()); //set the name of the value
    double tens[6]={0,0,0,0,0,0};
    for (int i=0;i<nelts;i++) {
     dodoF2C_GetElementValueTensor2s(i+1,cname,name.size(),tens);
     vals->InsertNextTupleValue(tens);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);

    delete[] cname;
    return 0;

}


int GetElementMaterialValueScalar(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts, std::string mat, std::string name, std::string newname)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());
    char *cmat = new char[mat.length() + 1]; // or
    std::strcpy(cmat, mat.c_str());



    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
    vals->SetName(newname.c_str()); //set the name of the value
    double scal=0.0;
    for (int i=0;i<nelts;i++) {
     dodoF2C_GetElementMaterialValueScalar(i+1,cname,name.size(),cmat,mat.size(),&scal);
     //std::cout << scal << std::endl << std::endl;
     vals->InsertNextValue(scal);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);

    delete[] cname;
    delete[] cmat;
    return 0;

}


int GetElementMaterialValueTensor1(vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid, int nelts, std::string mat, std::string name, std::string newname)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());
    char *cmat = new char[mat.length() + 1]; // or
    std::strcpy(cmat, mat.c_str());



    vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
    vals->SetNumberOfComponents(3); //we will have only 1 value associated with the triangle
    vals->SetName(newname.c_str()); //set the name of the value
    double tens[3]={0,0,0};
    for (int i=0;i<nelts;i++) {
     dodoF2C_GetElementMaterialValueTensor1(i+1,cname,name.size(),cmat,mat.size(),tens);
     //std::cout << scal << std::endl << std::endl;
     vals->InsertNextTupleValue(tens);
    }
    unstructuredGrid->GetCellData()->AddArray(vals);

    delete[] cname;
    delete[] cmat;
    return 0;

}


double GetScalar(int ivalue, std::string mat, std::string name)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());
    char *cmat = new char[mat.length() + 1]; // or
    std::strcpy(cmat, mat.c_str());

    double scal=0.0;
    
    dodoF2C_GetElementMaterialValueScalar(ivalue,cname,name.size(),cmat,mat.size(),&scal);

    delete[] cname;
    delete[] cmat;
    return scal;

}
tensor1 GetTensor1(int ivalue,  std::string mat, std::string name)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());
    char *cmat = new char[mat.length() + 1]; // or
    std::strcpy(cmat, mat.c_str());

    double tens[3]={0,0,0};

    dodoF2C_GetElementMaterialValueTensor1(ivalue,cname,name.size(),cmat,mat.size(),tens);
    tensor1 tensor1temp;
    tensor1temp=tens;
    delete[] cname;
    delete[] cmat;
    return tensor1temp;

}

tensor2 GetTensor2s(int ivalue, std::string name)
{

    char *cname = new char[name.length() + 1]; // or
    std::strcpy(cname, name.c_str());

    double tens[6]={0,0,0,0,0,0};

    dodoF2C_GetElementValueTensor2s(ivalue,cname,name.size(),tens);


    delete[] cname;
    tensor2 tensor2temp;
    tensor2temp=tens;
    return tensor2temp;

}


int main(int argc, char *argv[])
{

clock_t clstart,clend;
clstart=clock();
    
 int ierr = MPI_Init ( &argc, &argv );



  if(argc != 2)
  {
   std::cout << "Required arguments: OutputFilename" << std::endl;
   return EXIT_FAILURE;
    MPI_Finalize ( );
  }

  double ts0[6]={0,0,0,0,0,0};
   int nnodes,nelts,nframes;
  int nesets;
 
  int step,cyc,inc,it;
   double time, time_step, time_cycle;
   it=0;
   
   dodoF2C_Open(argv[1],strlen(argv[1]),&ierr);

   dodoF2C_GetNodesNumber(&nnodes);
   dodoF2C_GetElementsNumber(&nelts);
   dodoF2C_GetFramesNumber(&nframes);

      dodoF2C_GetEsetsNumber(&nesets);
      nesets=nesets-3;

   std::cout << nnodes << " " << nelts << " " << nframes << " "<<nesets << std::endl;

   
   
   
   double **tau_min_array=new double *[12];
   for(int i=0;i<12;i++)
   {    
       tau_min_array[i]= new double[nesets];
   }
   for(int i=0;i<12;i++)  
   {
       for(int j=0;j<nesets;j++)
       {
           tau_min_array[i][j]=100000.0;
       }
   }
   double **tau_max_array=new double *[12];
   for(int i=0;i<12;i++)
   {    
       tau_max_array[i]= new double[nesets];
   }
   for(int i=0;i<12;i++)  
   {
       for(int j=0;j<nesets;j++)
       {
           tau_max_array[i][j]=0.0;
       }
   }   
   double **tau_avg_array=new double *[12];
   for(int i=0;i<12;i++)
   {    
       tau_avg_array[i]= new double[nesets];
   }
   for(int i=0;i<12;i++)  
   {
       for(int j=0;j<nesets;j++)
       {
           tau_avg_array[i][j]=0.0;
       }
   }     
   
   
    double maxstress;
    maxstress=0.0; 
    tensor2 sigma_calculated,sigma_grain_calculated;
    tensor2 direction_system;
    tensor1 normal_direction_i,director_direction_i;    
    tensor1 tau_vector;
    double tau_value; 
    double tau_re;
    double stress;
    std::string normal_direction_name,director_direction_name,strI;
    double maxtau_re=0.0;
    int keystep,keycyc,keyinc,keyelement,keyeset,keyslipsystem;
    double esetinfo_volume,keytime_step,keytime_cycle,keytime;
    
   
    double hahahaha;
    dodoF2C_GetElementVolume(888, &hahahaha);
    std::cout<<"hahahaha"<<hahahaha<<std::endl;
    
    
   
   
   Elt *elt=new Elt[nelts];
   Grain *grain=new Grain[nesets];
   
   bool bb[nelts];
//   dodoF2C_GetElementsSet("grain1",6,bb);

 
   char *name = new char[20]; int namelength;double tempvolume;
   for (int n=0;n<nesets;n++)
   {
     dodoF2C_GetEsetName(n+1,name, 20);
        namelength=0;
        for(int ch=0;ch<20;ch++)
        {
            if (name[ch]!=' ') namelength+=1;
        }
        
      //  std::cout<<name<<namelength<<std::endl;
     dodoF2C_GetElementsSet(name,namelength,bb);
     //std::cout<<bb[1]<<std::endl;
     for(int i=0;i<nelts;i++)
     {
         if(bb[i]==255 && elt[i].f==false) 
         {
            // std::cout<<"good ";
             elt[i].gn=n+1;elt[i].f=true;
            dodoF2C_GetElementVolume(i+1, &tempvolume);
             elt[i].volume=tempvolume;
            grain[n].elts++;
            grain[n].n_elt.push_back(i);
             grain[n].volume=grain[n].volume+elt[i].volume;
        }
     }
    // std::cout<<"volume"<<grain[n].volume<<"      ";
     
   }

std::cout<<"first part goes"<<std::endl;
   
  
       double maxtauvalue=0.0;
   
for (int f=nframes/stabilized_cycle*(stabilized_cycle-1);f<nframes;f++)
   {
    dodoF2C_ReadFrame(f+1,&ierr);
    if (ierr!=0) break;
    dodoF2C_GetIncrement(&inc);
    dodoF2C_GetStep(&step);
    dodoF2C_GetCycle(&cyc);
    dodoF2C_GetTotalTime(&time);
    std::cout << time << std::endl;
    dodoF2C_GetStepTime(&time_step);
    dodoF2C_GetCycleTime(&time_cycle);

     std::cout<<"working"<<std::endl;
    
 
     
    for(int i=0;i<12;i++)
    {     
        std::stringstream ss;
        ss<<i+1;
        ss>>strI;
        ss.clear();
        normal_direction_name="ncs"+strI;
        director_direction_name="bcs"+strI;
    for(int k=0;k<nesets-1;k++)
    {
        sigma_grain_calculated=ts0;
        int temptemp;
        temptemp=grain[k].n_elt[0];
       
        normal_direction_i=GetTensor1(temptemp+1,"Mat4",normal_direction_name);
        director_direction_i=GetTensor1(temptemp+1,"Mat4",director_direction_name);
        direction_system=(kroneckerproduct(normal_direction_i,director_direction_i)+kroneckerproduct(director_direction_i,normal_direction_i))/2.0;
        
        
    for (int l=0;l<grain[k].elts;l++) 
        {
            int j;
            j=grain[k].n_elt[l];
            //std::cout<<"tt"<<std::endl;
            sigma_calculated=GetTensor2s(j+1,"sig");
            sigma_grain_calculated=sigma_calculated*(elt[j].volume/grain[k].volume)+sigma_grain_calculated;
            
        }
                tau_value=doubleinnerproduct(direction_system,sigma_grain_calculated);
                
                //if(tau_value>maxtauvalue) { maxtauvalue=tau_value ;std::cout<<"maxtauvalue"<<maxtauvalue<<std::endl; }
                
                if(tau_value>tau_max_array[i][k]) {tau_max_array[i][k]=tau_value;}
                
                if(tau_value<tau_min_array[i][k]) {tau_min_array[i][k]=tau_value;}
                
                tau_avg_array[i][k]=(tau_max_array[i][k]+tau_min_array[i][k])/2.0;
                
                //if (tau_avg_array[i][k]>50 ) std::cout<<i<<' '<<k<<' ' <<tau_avg_array[i][k]<<std::endl;
        }
        
    }
    
       std::cout<<"taumax"<<tau_max_array[0][1269]<<std::endl;
       std::cout<<"taumax"<<tau_min_array[0][1269]<<std::endl;
       std::cout<<"taumax"<<tau_avg_array[0][1269]<<std::endl;
   }
   std::cout<<"tau_avg read "<<std::endl;

   
    ierr=0;
    
    //std::cout<<tau_avg_array[11][4628]<<std::endl;
    //std::cout<<  normal_direction_i.a1<<normal_direction_i.a2<<normal_direction_i.a3  <<std::endl;
    //std::cout<<  director_direction_i.a1<<director_direction_i.a2<<director_direction_i.a3  <<std::endl;
    //std::cout<<  sigma_calculated.a11<<sigma_calculated.a22<<sigma_calculated.a33<<sigma_calculated.a21<<sigma_calculated.a13<<sigma_calculated.a23  <<std::endl;
    //std::cout<<  direction_system.a11<<direction_system.a22<<direction_system.a33<<direction_system.a12<<direction_system.a13<<direction_system.a23  <<std::endl;
    
    std::cout<<"second part begins"<<std::endl;
    
    
   for (int f=nframes/stabilized_cycle*(stabilized_cycle-1);f<nframes;f++)
   {
    dodoF2C_ReadFrame(f+1,&ierr);
    if (ierr!=0) break;
    dodoF2C_GetIncrement(&inc);
    dodoF2C_GetStep(&step);
    dodoF2C_GetCycle(&cyc);
    dodoF2C_GetTotalTime(&time);
    std::cout << time << std::endl;
    std::cout<<"working"<<std::endl;
    dodoF2C_GetStepTime(&time_step);
    dodoF2C_GetCycleTime(&time_cycle);

    for(int i=0;i<12;i++)
    {     
        std::stringstream ss;
        ss<<i+1;
        ss>>strI;
        ss.clear();
        normal_direction_name="ncs"+strI;
        director_direction_name="bcs"+strI;
    for(int k=0;k<nesets-1;k++)
    {
        sigma_grain_calculated=ts0;
        int temptemp;
        temptemp=grain[k].n_elt[0];
        normal_direction_i=GetTensor1(temptemp+1,"Mat4",normal_direction_name);
        director_direction_i=GetTensor1(temptemp+1,"Mat4",director_direction_name);
        direction_system=(kroneckerproduct(normal_direction_i,director_direction_i)+kroneckerproduct(director_direction_i,normal_direction_i))/2.0;
        
    for (int l=0;l<grain[k].elts;l++) 
        {
            int j;
            j=grain[k].n_elt[l];
            //std::cout<<"tt"<<std::endl;
            sigma_calculated=GetTensor2s(j+1,"sig");  
            sigma_grain_calculated=sigma_calculated*elt[j].volume/grain[k].volume+sigma_grain_calculated;
            
        }
        
         tau_value=doubleinnerproduct(direction_system,sigma_grain_calculated);
         
         tau_re=fabs(tau_value-tau_avg_array[i][k]);
                
         if ((tau_re-maxtau_re)>0.0) {maxtau_re=tau_re;  std::cout<<"maxtau-re"<<maxtau_re<<std::endl;      }
               
                stress=tau_re+alpha*sigmah(time)/3.0;
                
                
                
                if ((stress-maxstress)>0) 
                {
                    maxstress=stress;
                    keystep=step;
                    keycyc=cyc;
                    keyinc=inc;
                    keyeset=k+1;
                    keyslipsystem=i+1;
                    keytime=time;
                    keytime_step=time_step;
                    keytime_cycle=time_cycle;
                    
                    esetinfo_volume=grain[k].volume;
                    
                    
                }       
               
        }
        
    
    }
    
    

   }
   
   
   
   


   
   
   
  for(int i=0;i<12;i++)
  {
      delete[] tau_max_array[i];
      delete[] tau_min_array[i];
  }
  delete[] tau_max_array;
  delete[] tau_min_array; 
     
  for(int i=0;i<12;i++)
      delete[] tau_avg_array[i];
  delete[] tau_avg_array;

    std::ofstream filetxt;
     filetxt.open("result.txt");
     //filetxt<<"tau-resolu"<<"          "<<maxtau_re<<std::endl;
     filetxt<<"Dang Van stress"<<"          "<<maxstress<<std::endl;
     filetxt<<"keystep"<<"          "<<keystep<<std::endl;
     filetxt<<"keycyc"<<"           "<<keycyc<<std::endl;
     filetxt<<"keyinc"<<"           "<<keyinc<<std::endl;
     //filetxt<<"keyelement"<<"           "<<keyelement<<std::endl;
     filetxt<<"keyelementset"<<"           grain"<<keyeset<<std::endl;
     filetxt<<"keyslipsystem"<<"                "<<keyslipsystem<<std::endl;
     filetxt<<"keytime"<<"           "<<keytime<<std::endl;
     filetxt<<"keytime_step"<<"           "<<keytime_step<<std::endl;
     filetxt<<"keytime_cycle"<<"           "<<keytime_cycle<<std::endl;
     
    filetxt.close();
 
    std::cout<<std::endl;   
    std::cout<<"tau-resolu"<<"          "<<maxtau_re<<std::endl;
    std::cout<<"Dang Van stress"<<"          "<<maxstress<<std::endl;
    

  dodoF2C_Close();
  MPI_Finalize ( );
  
  clend=clock();
  
  std::cout<<"time used"<<(clend-clstart)/1000.0<<std::endl;
  
  

  

  
  return EXIT_SUCCESS;

/*
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(0, 0, 0);
  points->InsertNextPoint(1, 0, 0);
  points->InsertNextPoint(1, 1, 0);
  points->InsertNextPoint(0, 1, 1);

   vtkSmartPointer<vtkTetra> tetra =
    vtkSmartPointer<vtkTetra>::New();

  tetra->GetPointIds()->SetId(0, 0);
  tetra->GetPointIds()->SetId(1, 1);
  tetra->GetPointIds()->SetId(2, 2);
  tetra->GetPointIds()->SetId(3, 3);

  vtkSmartPointer<vtkCellArray> cellArray =
    vtkSmartPointer<vtkCellArray>::New();
  cellArray->InsertNextCell(tetra);


  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  unstructuredGrid->SetPoints(points);
  unstructuredGrid->SetCells(VTK_TETRA, cellArray);

  // Write file
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(unstructuredGrid);
#else
  writer->SetInputData(unstructuredGrid);
#endif
  
  vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
  vals->SetNumberOfComponents(1); //we will have only 1 value associated with the triangle
  vals->SetName("num"); //set the name of the value
  vals->InsertNextValue(1.0);
   //unstructuredGrid->Update(); 
  unstructuredGrid->GetCellData()->AddArray(vals);
  
  
  writer->Write();
*/
  
  /*
  // Read and display file for verification that it was written correclty
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkSmartPointer<vtkDataSetMapper> mapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetInputConnection(reader->GetOutputPort());

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(.3, .6, .3); // Background color green

  renderWindow->Render();
  renderWindowInteractor->Start();
  */

}
