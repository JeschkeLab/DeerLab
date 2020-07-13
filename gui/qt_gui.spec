# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['qt_gui.py'],
             pathex=['C:\\Windows\\System32\\downlevel', 'D:\\lufa\\projects\\DeerLab\\PyDeerLab\\PyDeerLab\\gui'],
             binaries=[],
             datas=[('qt_gui.ui','.')],
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
          [],
          exclude_binaries=True,
          name='qt_gui',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='qt_gui')
