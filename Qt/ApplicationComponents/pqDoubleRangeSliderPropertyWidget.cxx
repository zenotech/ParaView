/*=========================================================================

   Program: ParaView
   Module:    $RCSfile$

   Copyright (c) 2005,2006 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2.

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

========================================================================*/
#include "pqDoubleRangeSliderPropertyWidget.h"
#include "ui_pqDoubleRangeSliderPropertyWidget.h"

#include "pqPropertiesPanel.h"
#include "pqProxyWidget.h"
#include "pqDoubleRangeWidget.h"
#include "pqWidgetRangeDomain.h"
#include "vtkSMProperty.h"

#include <QGridLayout>

class pqDoubleRangeSliderPropertyWidget::pqInternals
{
public:
  Ui::DoubleRangeSliderPropertyWidget Ui;
};

//-----------------------------------------------------------------------------
pqDoubleRangeSliderPropertyWidget::pqDoubleRangeSliderPropertyWidget(
  vtkSMProxy* smProxy, vtkSMProperty* smProperty, QWidget* parentObject)
  : Superclass(smProxy, parentObject),
  Internals(new pqDoubleRangeSliderPropertyWidget::pqInternals())
{
  this->setShowLabel(false);
  this->setChangeAvailableAsChangeFinished(false);

  Ui::DoubleRangeSliderPropertyWidget &ui = this->Internals->Ui;
  ui.setupUi(this);

  ui.gridLayout->setMargin(pqPropertiesPanel::suggestedMargin());
  ui.gridLayout->setVerticalSpacing(pqPropertiesPanel::suggestedVerticalSpacing());
  ui.gridLayout->setHorizontalSpacing(pqPropertiesPanel::suggestedHorizontalSpacing());

  this->addPropertyLink(ui.ThresholdBetween_0, "value", SIGNAL(valueChanged(double)),
                        smProperty, 0);
  this->addPropertyLink(ui.ThresholdBetween_1, "value", SIGNAL(valueChanged(double)),
                        smProperty, 1);

  this->connect(ui.ThresholdBetween_0, SIGNAL(valueEdited(double)),
                   this, SLOT(lowerChanged(double)));
  this->connect(ui.ThresholdBetween_1, SIGNAL(valueEdited(double)),
                   this, SLOT(upperChanged(double)));

  /// pqWidgetRangeDomain ensures that whenever the domain changes, the slider's
  /// ranges are updated.
  new pqWidgetRangeDomain(ui.ThresholdBetween_0, "minimum", "maximum",
                          smProperty, 0);
  new pqWidgetRangeDomain(ui.ThresholdBetween_1, "minimum", "maximum",
                          smProperty, 1);
}

//-----------------------------------------------------------------------------
pqDoubleRangeSliderPropertyWidget::~pqDoubleRangeSliderPropertyWidget()
{
  delete this->Internals;
  this->Internals = NULL;
}

//-----------------------------------------------------------------------------
void pqDoubleRangeSliderPropertyWidget::lowerChanged(double val)
{
  // clamp the lower threshold if we need to
  if (this->Internals->Ui.ThresholdBetween_1->value() < val)
    {
    this->Internals->Ui.ThresholdBetween_1->setValue(val);
    }

  emit this->changeFinished();
}

//-----------------------------------------------------------------------------
void pqDoubleRangeSliderPropertyWidget::upperChanged(double val)
{
  // clamp the lower threshold if we need to
  if(this->Internals->Ui.ThresholdBetween_0->value() > val)
    {
    this->Internals->Ui.ThresholdBetween_0->setValue(val);
    }

  emit this->changeFinished();
}
