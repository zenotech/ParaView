// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
/**
 * @class   vtkNetworkImageSource
 * @brief   an image source that can read an image file on
 * one process and ensure that it's available on some other group of processes.
 *
 * vtkNetworkImageSource is a vtkImageAlgorithm that can read an image file on
 * on the client process and produce the output image data on a client and
 * render-server processes.
 */

#ifndef vtkNetworkImageSource_h
#define vtkNetworkImageSource_h

#include "vtkImageAlgorithm.h"
#include "vtkPVVTKExtensionsFiltersRenderingModule.h" //needed for exports

class vtkImageData;
class vtkClientServerStream;

class VTKPVVTKEXTENSIONSFILTERSRENDERING_EXPORT vtkNetworkImageSource : public vtkImageAlgorithm
{
public:
  static vtkNetworkImageSource* New();
  vtkTypeMacro(vtkNetworkImageSource, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  ///@{
  /**
   * Get/Set the filename.
   */
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  ///@}

  /**
   * Needs to be called to perform the actual image migration.
   */
  void UpdateImage();

protected:
  vtkNetworkImageSource();
  ~vtkNetworkImageSource() override;

  vtkTimeStamp UpdateImageTime;

  char* FileName;

  vtkImageData* Buffer;
  int ReadImageFromFile(const char* filename);
  int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;
  int RequestInformation(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

private:
  vtkNetworkImageSource(const vtkNetworkImageSource&) = delete;
  void operator=(const vtkNetworkImageSource&) = delete;
};

#endif
