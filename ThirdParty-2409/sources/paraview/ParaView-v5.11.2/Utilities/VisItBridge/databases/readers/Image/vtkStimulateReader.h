// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                           vtkStimulateReader.h                            //
// ************************************************************************* //

#ifndef VTK_STIMULATE_READER_H
#define VTK_STIMULATE_READER_H

#include <vtkImageReader2.h>


// ****************************************************************************
//  Class: vtkStimulateReader
//
//  Purpose:
//      Reads the "Stimulate" image format, which is what is generated by
//      the NDE folks in Engineering.
//
//  Programmer: Hank Childs
//  Creation:   March 17, 2005
//
//  Modifications:
//
//    Mark C. Miller, Fri Apr 23 23:32:27 PDT 2010
//    Added tdata and dataType.
// ****************************************************************************

class vtkStimulateReader : public vtkImageReader2
{
public:
  static vtkStimulateReader *New();
  vtkTypeMacro(vtkStimulateReader,vtkImageReader2);

  int CanReadFile(const char* fname) override;
  const char* GetFileExtensions() override
    {
      return ".sdt .spr .SDT .SPR";
    }

  const char* GetDescriptiveName() override
    {
      return "Stimulate";
    }
  
  void GetDimensions(int &x, int &y) { x = dims[0]; y = dims[1]; };
  void GetOrigin(float &x, float &y) { x = origin[0]; y = origin[1]; };
  void GetStep(float &x, float &y) { x = step[0]; y = step[1]; };


protected:
  enum dtype {UCHAR, SHORT, INT, FLOAT};
  int RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector) override;
  void ExecuteDataWithInformation(vtkDataObject *, vtkInformation* outInfo) override;
  int OpenFile(void) override;
  vtkStimulateReader();
  ~vtkStimulateReader();
private:
  bool     haveReadSPRFile;
  bool     validSPRFile;
  int      dims[2];
  float    origin[2];
  float    step[2];
  dtype    dataType;

  vtkStimulateReader(const vtkStimulateReader&);  // Not implemented.
  void operator=(const vtkStimulateReader&);  // Not implemented.

  bool GetFilenames(const char *, char *spr_name, char *sdt_name);
  bool ReadSPRFile(const char *);
};

#endif

