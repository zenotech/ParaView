// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-FileCopyrightText: Copyright (c) Sandia Corporation
// SPDX-License-Identifier: BSD-3-Clause
#include "pqLiveSourceBehavior.h"

#include "pqApplicationCore.h"
#include "pqPipelineSource.h"
#include "pqServer.h"
#include "pqServerManagerModel.h"
#include "pqTimer.h"
#include "pqView.h"
#include "vtkClientServerStream.h"
#include "vtkCommand.h"
#include "vtkPVXMLElement.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSMSession.h"
#include "vtkSMSourceProxy.h"
#include "vtkSMViewProxy.h"

#include <QPointer>
#include <algorithm>
#include <cassert>
#include <limits>
#include <utility>
#include <vector>

bool pqLiveSourceBehavior::PauseLiveUpdates = false;

//-----------------------------------------------------------------------------
class pqLiveSourceBehavior::pqInternals
{
  void updateSources()
  {
    const auto old_size = this->LiveSources.size();
    auto iter = this->LiveSources.begin();
    while (iter != this->LiveSources.end())
    {
      if (iter->first == nullptr)
      {
        iter = this->LiveSources.erase(iter);
      }
      else
      {
        ++iter;
      }
    }

    if (old_size != this->LiveSources.size())
    {
      // update interval
      int interval = std::numeric_limits<int>::max();
      for (const auto& pair : this->LiveSources)
      {
        interval = std::min(interval, pair.second);
      }
      this->Timer.setInterval(interval);
    }
  }

  static const int DEFAULT_INTERVAL = 100;

public:
  pqTimer Timer;
  std::vector<std::pair<QPointer<pqPipelineSource>, int>> LiveSources;

  pqInternals()
  {
    this->Timer.setInterval(std::numeric_limits<int>::max());
    this->Timer.setSingleShot(true);
  }

  void addSource(pqPipelineSource* src, vtkPVXMLElement* liveHints)
  {
    assert(liveHints != nullptr);
    int interval = 0;
    if (!liveHints->GetScalarAttribute("interval", &interval) || interval <= 0)
    {
      interval = DEFAULT_INTERVAL;
    }

    this->LiveSources.push_back(std::make_pair(QPointer<pqPipelineSource>(src), interval));
    this->Timer.setInterval(std::min(this->Timer.interval(), interval));
    this->Timer.start();
  }

  void pause() { this->Timer.stop(); }

  void resume()
  {
    this->updateSources();
    if (!this->LiveSources.empty())
    {
      this->Timer.start();
    }
  }

  void update()
  {
    this->Timer.stop();
    this->updateSources();
    if (this->LiveSources.empty())
    {
      return;
    }

    if (this->tryAgainLater())
    {
      this->Timer.start();
      return;
    }

    // iterate over all sources and update those that need updating.
    if (!pqLiveSourceBehavior::isPaused())
    {
      for (const auto& pair : this->LiveSources)
      {
        pqPipelineSource* src = pair.first;
        auto proxy = vtkSMSourceProxy::SafeDownCast(src->getProxy());
        auto session = proxy->GetSession();

        vtkClientServerStream stream;
        stream << vtkClientServerStream::Invoke << VTKOBJECT(proxy) << "GetNeedsUpdate"
               << vtkClientServerStream::End;
        session->ExecuteStream(vtkPVSession::DATA_SERVER_ROOT, stream, /*ignore errors*/ true);

        vtkClientServerStream result = session->GetLastResult(vtkPVSession::DATA_SERVER_ROOT);
        bool needs_update = false;
        if (result.GetNumberOfMessages() == 1 && result.GetNumberOfArguments(0) == 1 &&
          result.GetArgument(0, 0, &needs_update) && needs_update)
        {
          proxy->MarkModified(proxy);
          src->renderAllViews();
        }
      }
    }

    this->Timer.start();
  }

private:
  bool tryAgainLater()
  {
    pqServerManagerModel* smmodel = pqApplicationCore::instance()->getServerManagerModel();
    auto servers = smmodel->findItems<pqServer*>();
    for (auto server : servers)
    {
      if (server->isProcessingPending() || server->session()->GetPendingProgress())
      {
        return true;
      }
    }
    return false;
  }
};

//-----------------------------------------------------------------------------
pqLiveSourceBehavior::pqLiveSourceBehavior(QObject* parentObject)
  : Superclass(parentObject)
  , Internals(new pqLiveSourceBehavior::pqInternals())
{
  pqServerManagerModel* smmodel = pqApplicationCore::instance()->getServerManagerModel();
  this->connect(smmodel, SIGNAL(viewAdded(pqView*)), SLOT(viewAdded(pqView*)));
  this->connect(
    smmodel, SIGNAL(sourceAdded(pqPipelineSource*)), SLOT(sourceAdded(pqPipelineSource*)));

  this->connect(&this->Internals->Timer, SIGNAL(timeout()), SLOT(timeout()));

  Q_FOREACH (pqView* view, smmodel->findItems<pqView*>())
  {
    this->viewAdded(view);
  }

  Q_FOREACH (pqPipelineSource* src, smmodel->findItems<pqPipelineSource*>())
  {
    this->sourceAdded(src);
  }
}

//-----------------------------------------------------------------------------
pqLiveSourceBehavior::~pqLiveSourceBehavior() = default;

//-----------------------------------------------------------------------------
void pqLiveSourceBehavior::pause()
{
  pqLiveSourceBehavior::PauseLiveUpdates = true;
}

//-----------------------------------------------------------------------------
void pqLiveSourceBehavior::resume()
{
  pqLiveSourceBehavior::PauseLiveUpdates = false;
}

//-----------------------------------------------------------------------------
void pqLiveSourceBehavior::viewAdded(pqView* view)
{
  if (auto viewProxy = vtkSMViewProxy::SafeDownCast(view->getProxy()))
  {
    if (auto iren = viewProxy->GetInteractor())
    {
      iren->AddObserver(
        vtkCommand::StartInteractionEvent, this, &pqLiveSourceBehavior::startInteractionEvent);
      iren->AddObserver(
        vtkCommand::EndInteractionEvent, this, &pqLiveSourceBehavior::endInteractionEvent);
    }
  }
}

//-----------------------------------------------------------------------------
void pqLiveSourceBehavior::sourceAdded(pqPipelineSource* src)
{
  vtkSMProxy* proxy = src->getProxy();
  if (auto hints = proxy->GetHints())
  {
    if (auto liveHints = hints->FindNestedElementByName("LiveSource"))
    {
      this->Internals->addSource(src, liveHints);
    }
  }
}

//-----------------------------------------------------------------------------
void pqLiveSourceBehavior::startInteractionEvent()
{
  this->Internals->pause();
}

//-----------------------------------------------------------------------------
void pqLiveSourceBehavior::endInteractionEvent()
{
  this->Internals->resume();
}

//-----------------------------------------------------------------------------
void pqLiveSourceBehavior::timeout()
{
  this->Internals->update();
}
