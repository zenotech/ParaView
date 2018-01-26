/*=========================================================================

  Program:   ParaView
  Module:    vtkChartRepresentation.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkChartRepresentation
 *
 * vtkChartRepresentation is the base representation for charting
 * representations. Currently, ParaView's charting views are client-side only
 * views that render only on the client side. That being the case, when running
 * in client-server mode or in parallel, the data-delivery mode is fixed. Hence,
 * unlike representations for 3D views, this representation delivers the data in
 * RequestData() itself. This makes it possible for client code to call
 * UpdatePipeline() on the representation proxy and then access the delivered
 * vtkTable on the client.
*/

#ifndef vtkChartRepresentation_h
#define vtkChartRepresentation_h

#include "vtkPVDataRepresentation.h"
#include "vtkSmartPointer.h" // needed for vtkSmartPointer
#include "vtkStdString.h"    //  needed for vtkStdString.
#include "vtkWeakPointer.h"  // needed for vtkWeakPointer

#include <map> // needed for map
#include <set> //needed for ivars

class vtkChartSelectionRepresentation;
class vtkCSVExporter;
class vtkMultiBlockDataSet;
class vtkPVCacheKeeper;
class vtkPVContextView;
class vtkSelectionDeliveryFilter;
class vtkTable;

class VTKPVCLIENTSERVERCORERENDERING_EXPORT vtkChartRepresentation : public vtkPVDataRepresentation
{
public:
  static vtkChartRepresentation* New();
  vtkTypeMacro(vtkChartRepresentation, vtkPVDataRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  /**
   * These must only be set during initialization before adding the
   * representation to any views or calling Update().
   */
  void SetSelectionRepresentation(vtkChartSelectionRepresentation*);

  /**
   * Set visibility of the representation.
   */
  virtual void SetVisibility(bool visible) VTK_OVERRIDE;

  /**
   * This needs to be called on all instances of vtkGeometryRepresentation when
   * the input is modified. This is essential since the geometry filter does not
   * have any real-input on the client side which messes with the Update
   * requests.
   */
  virtual void MarkModified() VTK_OVERRIDE;

  // *************************************************************************

  //@{
  /**
   * Set the field association for arrays to use. When changed, this will call
   * MarkModified().
   */
  void SetFieldAssociation(int);
  vtkGetMacro(FieldAssociation, int);
  //@}

  // methods to control block selection.
  // When changed, this will call MarkModified().
  void SetCompositeDataSetIndex(unsigned int); // only used for single block selection
  void AddCompositeDataSetIndex(unsigned int);
  void ResetCompositeDataSetIndices();

  /**
   * Override because of internal selection representations that need to be
   * initialized as well.
   */
  virtual unsigned int Initialize(
    unsigned int minIdAvailable, unsigned int maxIdAvailable) VTK_OVERRIDE;

  /**
   * vtkAlgorithm::ProcessRequest() equivalent for rendering passes. This is
   * typically called by the vtkView to request meta-data from the
   * representations or ask them to perform certain tasks e.g.
   * PrepareForRendering.
   * Overridden to handle REQUEST_RENDER() to call PrepareForRendering.
   */
  virtual int ProcessViewRequest(vtkInformationRequestKey* request_type, vtkInformation* inInfo,
    vtkInformation* outInfo) VTK_OVERRIDE;

  /**
   * Method to provide the default name given the name of a table and a column
   * in that table. When overriding this method, ensure that the implementation
   * does not depend on the instance being initialized or have valid data. This
   * method is called on the client-side when filling up the domains for the
   * properties.
   */
  virtual vtkStdString GetDefaultSeriesLabel(
    const vtkStdString& tableName, const vtkStdString& columnName);

  //@{
  /**
   * Flatten the table, i.e. split any multicomponent columns into separate
   * components.
   */
  vtkSetMacro(FlattenTable, int);
  vtkGetMacro(FlattenTable, int);
  //@}

  /**
   * This method is called on the client-side by the vtkPVContextView whenever a
   * new selection is made on all the visible representations in that view.
   * The goal of this method is allow the representations to transform the
   * selection created in the view (which is an id-based selection based on the
   * vtkTable that is fed into the vtkChart) to an appropriate selection based
   * on the data going into the representation.
   * Return false if the selection is not applicable to this representation or
   * the conversion cannot be made.
   * Default implementation simply ensures that the FieldType on the
   * selection nodes is set to match up with the FieldAssociation on the
   * representation.
   */
  virtual bool MapSelectionToInput(vtkSelection* sel);

  /**
   * This is the inverse of MapSelectionToInput(). In this case, we are
   * converting the selection defined on the input for the representation to a
   * selection that corresponds to elements being rendered in the view.
   * The default implementation checks removes vtkSelectionNode items in sel
   * that don't have the FieldType that matches this->FieldAssociation.
   * Similar to MapSelectionToInput(), this method is expected to transform the
   * sel in place and return false is the selection is not applicable to this
   * representation or the conversion cannot be made.
   */
  virtual bool MapSelectionToView(vtkSelection* sel);

  /**
   * Called by vtkPVContextView::Export() to export the representation's data to
   * a CSV file. Return false on failure which will call the exporting process
   * to abort and raise an error. Default implementation simply returns false.
   */
  virtual bool Export(vtkCSVExporter* vtkNotUsed(exporter)) { return false; }

protected:
  vtkChartRepresentation();
  ~vtkChartRepresentation();

  /**
   * Fill input port information.
   */
  virtual int FillInputPortInformation(int port, vtkInformation* info) VTK_OVERRIDE;

  /**
   * This method is called before actual render if this->MTime was modified
   * since the last time this method was called. Subclasses should override to
   * update "appearance" related changes that don't affect data.
   * When this method is called, you're assured that this->ContextView is
   * valid.
   * Note that this method will not be called if this->GetVisibility()
   * returns false, this subclasses should also override SetVisibility() to
   * hide "actors" and such.
   */
  virtual void PrepareForRendering() {}

  /**
   * Subclasses should override this to connect inputs to the internal pipeline
   * as necessary. Since most representations are "meta-filters" (i.e. filters
   * containing other filters), you should create shallow copies of your input
   * before connecting to the internal pipeline. The convenience method
   * GetInternalOutputPort will create a cached shallow copy of a specified
   * input for you. The related helper functions GetInternalAnnotationOutputPort,
   * GetInternalSelectionOutputPort should be used to obtain a selection or
   * annotation port whose selections are localized for a particular input data object.
   */
  virtual int RequestData(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*) VTK_OVERRIDE;

  /**
   * Adds the representation to the view.  This is called from
   * vtkView::AddRepresentation().  Subclasses should override this method.
   * Returns true if the addition succeeds.
   */
  virtual bool AddToView(vtkView* view) VTK_OVERRIDE;

  /**
   * Removes the representation to the view.  This is called from
   * vtkView::RemoveRepresentation().  Subclasses should override this method.
   * Returns true if the removal succeeds.
   */
  virtual bool RemoveFromView(vtkView* view) VTK_OVERRIDE;

  /**
   * Overridden to check with the vtkPVCacheKeeper to see if the key is cached.
   */
  virtual bool IsCached(double cache_key) VTK_OVERRIDE;

  /**
   * Convenience method to get the first vtkTable from LocalOutput, if any.
   */
  vtkTable* GetLocalOutput();

  /**
   * Method to be overrided to transform input data to a vtkTable.
   * The default implementation just returns the data object provided in parameter.
   */
  virtual vtkDataObject* TransformInputData(
    vtkInformationVector** inputVector, vtkDataObject* data);

  /**
   * Method to be overridden to apply an operation of the table after it is
   * gathered to the first rank for rendering the chart.  This allows subclasses
   * to operate on the final table.  The default implementation just returns the
   * input.
   */
  virtual vtkSmartPointer<vtkDataObject> TransformTable(vtkSmartPointer<vtkDataObject> table);

  typedef std::map<std::string, vtkSmartPointer<vtkTable> > MapOfTables;
  /**
   * Convenience method to get all vtkTable instances with their associated
   * names.
   */
  bool GetLocalOutput(MapOfTables& tables);

  int FieldAssociation;
  vtkPVCacheKeeper* CacheKeeper;
  vtkWeakPointer<vtkPVContextView> ContextView;
  bool EnableServerSideRendering;
  int FlattenTable;

  vtkSmartPointer<vtkMultiBlockDataSet> LocalOutput;
  std::set<unsigned int> CompositeIndices; // the selected blocks

  vtkWeakPointer<vtkChartSelectionRepresentation> SelectionRepresentation;

private:
  vtkChartRepresentation(const vtkChartRepresentation&) VTK_DELETE_FUNCTION;
  void operator=(const vtkChartRepresentation&) VTK_DELETE_FUNCTION;

  vtkTimeStamp PrepareForRenderingTime;
  vtkSmartPointer<vtkChartSelectionRepresentation> DummyRepresentation;
};

#endif
