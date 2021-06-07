/*=========================================================================

  Program:   ParaView
  Module:    vtkTransferFunctionChartHistogram2D.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef vtkTransferFunctionChartHistogram2D_h
#define vtkTransferFunctionChartHistogram2D_h

#include "vtkChartHistogram2D.h"

#include "vtkCommand.h"             // needed for vtkCommand::UserEvent
#include "vtkRemotingViewsModule.h" // needed for export macro
#include "vtkWeakPointer.h"         // needed for vtkWeakPointer

// Forward declarations
class vtkContextMouseEvent;
class vtkImageData;
class vtkTransferFunctionBoxItem;

class VTKREMOTINGVIEWS_EXPORT vtkTransferFunctionChartHistogram2D : public vtkChartHistogram2D
{
public:
  static vtkTransferFunctionChartHistogram2D* New();
  vtkTypeMacro(vtkTransferFunctionChartHistogram2D, vtkChartHistogram2D);

  // Events fires by this class (and subclasses).
  // \li TransferFunctionModified is fired when the 2D transfer function is modified.
  enum
  {
    TransferFunctionModified = vtkCommand::UserEvent + 1000,
  };

  /**
   * Get whether the chart has been initialized with a histogram.
   */
  bool IsInitialized();

  /**
   * Add a new box item to the chart
   */
  vtkSmartPointer<vtkTransferFunctionBoxItem> AddNewBox();

  /**
   * Override to add new box item to the chart when double clicked.
   */
  bool MouseDoubleClickEvent(const vtkContextMouseEvent& mouse) override;

  /**
   * Set the input image data for the 2D histogram.
   */
  void SetInputData(vtkImageData*, vtkIdType z = 0) override;

  ///@{
  /**
   * Set/Get the image data that will be populated with the 2D transfer function.
   */
  virtual void SetTransferFunction2D(vtkImageData* transfer2D);
  virtual vtkImageData* GetTransferFunction2D();
  ///@}

protected:
  vtkTransferFunctionChartHistogram2D() {}
  ~vtkTransferFunctionChartHistogram2D() override = default;

  /**
   * Update individual item bounds based on the chart range.
   */
  void UpdateItemsBounds(
    const double xMin, const double xMax, const double yMin, const double yMax);

  /**
   * Generate the 2D transfer function from the box items.
   */
  virtual void GenerateTransfer2D();

  /**
   * Callback listening to the SelectionChangedEvent of box items to indicate that the 2D transfer
   * function was modified.
   */
  void OnTransferFunctionBoxItemModified(vtkObject* caller, unsigned long eid, void* callData);

  // Member variables;
  vtkWeakPointer<vtkImageData> TransferFunction2D;

private:
  vtkTransferFunctionChartHistogram2D(const vtkTransferFunctionChartHistogram2D&);
  void operator=(const vtkTransferFunctionChartHistogram2D&);
};

#endif // vtkTransferFunctionChartHistogram2D_h
