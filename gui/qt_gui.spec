# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['qt_gui.py'],
             pathex=['C:\\Windows\\System32\\downlevel', 'D:\\lufa\\projects\\DeerLab\\PyDeerLab\\PyDeerLab\\gui'],
             binaries=[],
             datas=[('qt_gui.ui', '.')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='qt_gui',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=['vcruntime140.dll', 'ucrtbase.dll', 'Qt5Core.dll', 'Qt5DBus.dll', 'Qt5Gui.dll', 'Qt5Network.dll', 'Qt5OpenGL.dll', 'Qt5Qml.dll', 'Qt5QmlModels.dll', 'Qt5Quick.dll', 'Qt5Svg.dll', 'Qt5Test.dll', 'Qt5WebSockets.dll', 'Qt5Widgets.dll'],
          runtime_tmpdir=None,
          console=True )
