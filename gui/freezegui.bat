
call pyinstaller --onefile --noconfirm --clean^
                 --add-data qt_gui.ui;.^
                 --paths C:\Windows\System32\downlevel^
                 --upx-dir=D:\lufa\code\upx-3.96-win64^
                 --upx-exclude=vcruntime140.dll^
                 --upx-exclude=ucrtbase.dll^
                 --upx-exclude=Qt5Core.dll^
                 --upx-exclude=Qt5DBus.dll^
                 --upx-exclude=Qt5Gui.dll^
                 --upx-exclude=Qt5Network.dll^
                 --upx-exclude=Qt5OpenGL.dll^
                 --upx-exclude=Qt5Qml.dll^
                 --upx-exclude=Qt5QmlModels.dll^
                 --upx-exclude=Qt5Quick.dll^
                 --upx-exclude=Qt5Svg.dll^
                 --upx-exclude=Qt5Test.dll^
                 --upx-exclude=Qt5WebSockets.dll^
                 --upx-exclude=Qt5Widgets.dll^
                  qt_gui.py



:: --add-data qt_gui.ui;. (Windows)  --add-data qt_gui.ui:. (Unix)
