/*=========================================================================

   Program: ParaView
   Module:  pqOutputWidget.h

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
#ifndef pqOutputWidget_h
#define pqOutputWidget_h

#include <QWidget>

#include "pqCoreModule.h" // for exports
#include <QScopedPointer> // for QScopedPointer.

class pqOutputWidget;

namespace OutputWidgetInternals
{
class OutputWindow;
};

Q_DECLARE_METATYPE(QtMsgType)
/**
 * @class MessageHandler
 * @brief widget to show Qt and VTK messages.
 *
 * This utility class allows us to install a Qt message handler. We are using a
 * signal/slot to ensure that the messages are written on the main thread.
 *
 * This utility class also helps us support VTK messages coming from different
 * threads.
 *
 * N.B. This is really a internal class however has to be exported so it can
 * be used inside a message handler function inside Qt.
 */
class PQCORE_EXPORT MessageHandler : public QObject
{
  Q_OBJECT
public:
  MessageHandler(QObject* parent = nullptr);

  /** This function can be used to install a message handler for Qt that forwards
   * all messages generated by Qt to the vtkOutputWindow instance obtained by
   * vtkOutputWindow::GetInstance(). Thus, this is quick way to forward all Qt
   * messages to be forwarded to VTK for output using the same mechanisms setup
   * for handling VTK messages.
   */
  static void install(pqOutputWidget* widget);

  /**
   * This is the function installed with Qt.
   */
  static void handler(QtMsgType type, const QMessageLogContext&, const QString& msg);

  /**
   * This function is used for VTK messages that don't have QMessageLogContext.
   */
  static void handlerVTK(QtMsgType type, const QString& msg);

Q_SIGNALS:
  void message(QtMsgType type, const QString& msg);
  void showMessage(QString msg, QtMsgType type);

private:
  static MessageHandler* instance();

private Q_SLOTS:
  void displayMessage(QtMsgType type, const QString& msg);
};

/**
 * @class pqOutputWidget
 * @brief widget to show Qt and VTK messages.
 *
 * pqOutputWidget is QWidget subclass that can capture and show messages
 * generated by Qt as well as VTK. On instantiation, it installs handlers to
 * handle Qt and VTK messages.
 *
 */
class PQCORE_EXPORT pqOutputWidget : public QWidget
{
  Q_OBJECT
  typedef QWidget Superclass;
  Q_PROPERTY(QString settingsKey READ settingsKey WRITE setSettingsKey);

public:
  pqOutputWidget(QWidget* parent = nullptr, Qt::WindowFlags f = Qt::WindowFlags());
  ~pqOutputWidget() override;

  /**
   * Add substrings to match with message to determine whether they should be
   * suppressed.
   */
  void suppress(const QStringList& substrs);

  //@{
  /**
  * Set the settings key that will be used to restore/save the state for
  * the "Show full messages" button. By default "pqOutputWidget" will be used.
  */
  void setSettingsKey(const QString& key);
  const QString& settingsKey() const;
  //@}

  //@{
  /**
   * Set/get the font size in points for the output console text.
   */
  void setFontSize(int fontSize);
  //@}

public Q_SLOTS:
  /**
   * Display a message in the widget. There's generally no need to use this
   * since one can have the same effect by triggering the message through Qt or
   * VTK message notification mechanisms.
   *
   * @returns true if the message was displayed, otherwise false.
   */
  bool displayMessage(const QString& message, QtMsgType type = QtInfoMsg);

  /**
   * Show full messages instead of grouped messages.
   */
  void showFullMessages(bool val);

  /**
   * Copy output text to the clipboard
   */
  void copyToClipboard();

  /**
   * Saves output text to a file.
   */
  void saveToFile();

  /**
   * Clears the console.
   */
  void clear();

Q_SIGNALS:
  /**
   * Signal fired any time a message is displayed.
   */
  void messageDisplayed(const QString& message, int type);

protected:
  /**
   * Returns `true` if the message must be suppressed/ignored.
   */
  virtual bool suppress(const QString& message, QtMsgType type);

  /**
   * Extract a summary string from the message and returns that.
   */
  virtual QString extractSummary(const QString& message, QtMsgType type);

private:
  Q_DISABLE_COPY(pqOutputWidget)

  class pqInternals;
  QScopedPointer<pqInternals> Internals;
  friend class OutputWidgetInternals::OutputWindow;
};

#endif
