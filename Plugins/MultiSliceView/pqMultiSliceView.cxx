/*=========================================================================

  Program:   Visualization Toolkit
  Module:    pqMultiSliceView.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "pqMultiSliceView.h"

#include <QtCore>
#include <QtGui>

#include "QVTKWidget.h"

#include "pqMultiSliceAxisWidget.h"
#include "pqRepresentation.h"

#include "vtkAxis.h"
#include "vtkChartXY.h"
#include "vtkDataRepresentation.h"
#include "vtkMath.h"
#include "vtkPVDataInformation.h"
#include "vtkPVRenderView.h"
#include "vtkPlot.h"
#include "vtkSMPropertyHelper.h"
#include "vtkSMRenderViewProxy.h"
#include "vtkSMRepresentationProxy.h"
#include "vtkView.h"

#define MULTI_SLICE_AXIS_THIKNESS 60
#define MULTI_SLICE_AXIS_ACTIVE_SIZE 20
#define MULTI_SLICE_AXIS_EDGE_MARGIN 10

//-----------------------------------------------------------------------------
pqMultiSliceView::pqMultiSliceView(
    const QString& viewType, const QString& group, const QString& name,
    vtkSMViewProxy* viewProxy, pqServer* server, QObject* p)
  : pqRenderView(viewType, group, name, viewProxy, server, p)
{
  QObject::connect(this, SIGNAL(representationAdded(pqRepresentation*)),
                   this, SLOT(updateAxisBounds()));
  QObject::connect(this, SIGNAL(representationRemoved(pqRepresentation*)),
                   this, SLOT(updateAxisBounds()));
  QObject::connect(this, SIGNAL(representationVisibilityChanged(pqRepresentation*,bool)),
                   this, SLOT(updateAxisBounds()));

  // Make sure all the representations share the same slice values
  QObject::connect(this, SIGNAL(representationAdded(pqRepresentation*)),
                   this, SLOT(updateSlices()));
}

//-----------------------------------------------------------------------------
pqMultiSliceView::~pqMultiSliceView()
{
}
//-----------------------------------------------------------------------------
QWidget* pqMultiSliceView::createWidget()
{
  // Get the internal widget that we want to decorate
  this->InternalWidget = qobject_cast<QVTKWidget*>(this->pqRenderView::createWidget());

  // Build the widget hierarchy
  QWidget* container = new QWidget();
  container->setStyleSheet("background-color: white");
  container->setAutoFillBackground(true);

  QGridLayout* gridLayout = new QGridLayout(container);
  this->InternalWidget->setParent(container);

  // Init top axis
  this->AxisX = new pqMultiSliceAxisWidget(container);
  this->AxisX->setAxisType(vtkAxis::LEFT);
  this->AxisX->setRange(-10,10);
  this->AxisX->setTitle("X");
  this->AxisX->SetEdgeMargin(MULTI_SLICE_AXIS_EDGE_MARGIN);
  this->AxisX->SetActiveSize(MULTI_SLICE_AXIS_ACTIVE_SIZE);
  this->AxisX->setFixedWidth(MULTI_SLICE_AXIS_THIKNESS);
  this->AxisX->renderView();

  this->AxisY = new pqMultiSliceAxisWidget(container);
  this->AxisY->setAxisType(vtkAxis::TOP);
  this->AxisY->setRange(-10,10);
  this->AxisY->setTitle("Y");
  this->AxisY->SetEdgeMargin(MULTI_SLICE_AXIS_EDGE_MARGIN);
  this->AxisY->SetActiveSize(MULTI_SLICE_AXIS_ACTIVE_SIZE);
  this->AxisY->setFixedHeight(MULTI_SLICE_AXIS_THIKNESS - 4);
  this->AxisY->renderView();

  this->AxisZ = new pqMultiSliceAxisWidget(container);
  this->AxisZ->setAxisType(vtkAxis::RIGHT);
  this->AxisZ->setRange(-10,10);
  this->AxisZ->setTitle("Z");
  this->AxisZ->SetEdgeMargin(MULTI_SLICE_AXIS_EDGE_MARGIN);
  this->AxisZ->SetActiveSize(MULTI_SLICE_AXIS_ACTIVE_SIZE);
  this->AxisZ->setFixedWidth(MULTI_SLICE_AXIS_THIKNESS);
  this->AxisZ->renderView();

  gridLayout->addWidget(this->AxisY, 0, 1);  // TOP
  gridLayout->addWidget(this->AxisX, 1, 0); // LEFT
  gridLayout->addWidget(this->AxisZ, 1, 2); // RIGHT
  gridLayout->addWidget(this->InternalWidget,1,1);
  gridLayout->setContentsMargins(0,0,0,0);
  gridLayout->setSpacing(0);

  // Properly do the binding between the proxy and the 3D widget
  vtkSMRenderViewProxy* renModule = this->getRenderViewProxy();
  if (this->InternalWidget && renModule)
    {
    this->InternalWidget->SetRenderWindow(renModule->GetRenderWindow());
    }

  // Make sure we are aware when the slice positions changes
  QObject::connect(this->AxisX, SIGNAL(modelUpdated()),
                   this, SLOT(updateSlices()));
  QObject::connect(this->AxisY, SIGNAL(modelUpdated()),
                   this, SLOT(updateSlices()));
  QObject::connect(this->AxisZ, SIGNAL(modelUpdated()),
                   this, SLOT(updateSlices()));

  return container;
}
//-----------------------------------------------------------------------------
void pqMultiSliceView::updateAxisBounds(double bounds[6])
{
  this->AxisX->setRange(bounds[0], bounds[1]);
  this->AxisY->setRange(bounds[2], bounds[3]);
  this->AxisZ->setRange(bounds[4], bounds[5]);

  // Make sure we render the new range
  this->AxisX->renderView();
  this->AxisY->renderView();
  this->AxisZ->renderView();
}

//-----------------------------------------------------------------------------
void pqMultiSliceView::updateAxisBounds()
{
  double bounds[6] = { VTK_DOUBLE_MAX, VTK_DOUBLE_MIN,
                       VTK_DOUBLE_MAX, VTK_DOUBLE_MIN,
                       VTK_DOUBLE_MAX, VTK_DOUBLE_MIN};
  foreach(pqRepresentation* rep, this->getRepresentations())
    {
    if(!rep->isVisible() || rep->isWidgetType())
      {
      continue;
      }
    vtkSMRepresentationProxy* smRep =
        vtkSMRepresentationProxy::SafeDownCast(rep->getProxy());
    smRep->UpdatePipeline();
    double* repBounds = smRep->GetRepresentedDataInformation()->GetBounds();
    for(int i=0;i<3;i++)
      {
      int index = 2*i;
      bounds[index] = (bounds[index] < repBounds[index]) ? bounds[index] : repBounds[index];
      index++;
      bounds[index] = (bounds[index] > repBounds[index]) ? bounds[index] : repBounds[index];
      }
    }

  if(bounds[0] != VTK_DOUBLE_MAX)
    {
    this->updateAxisBounds(bounds);
    }
}
//-----------------------------------------------------------------------------
void pqMultiSliceView::updateSlices()
{
  int nbSliceX = 0;
  const double* sliceX = this->AxisX->getVisibleSlices(nbSliceX);
  int nbSliceY = 0;
  const double* sliceY = this->AxisY->getVisibleSlices(nbSliceY);
  int nbSliceZ = 0;
  const double* sliceZ = this->AxisZ->getVisibleSlices(nbSliceZ);

  foreach(pqRepresentation* rep, this->getRepresentations())
    {
    if( !rep->isVisible() || rep->isWidgetType() ||
        !rep->getProxy()->GetProperty("XSlicesValues"))
      {
      continue;
      }
    vtkSMRepresentationProxy* smRep =
        vtkSMRepresentationProxy::SafeDownCast(rep->getProxy());

    vtkSMPropertyHelper(smRep, "XSlicesValues").Set(sliceX, nbSliceX);
    vtkSMPropertyHelper(smRep, "YSlicesValues").Set(sliceY, nbSliceY);
    vtkSMPropertyHelper(smRep, "ZSlicesValues").Set(sliceZ, nbSliceZ);

    smRep->MarkModified(NULL);
    smRep->MarkDirty(NULL);
    smRep->UpdateVTKObjects();
    }
  this->render();
}
