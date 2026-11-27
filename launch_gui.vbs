Option Explicit

Dim files, shell, project, interpreter, command, argument
Set files = CreateObject("Scripting.FileSystemObject")
Set shell = CreateObject("WScript.Shell")
project = files.GetParentFolderName(WScript.ScriptFullName)
interpreter = files.BuildPath(project, ".venv\Scripts\pythonw.exe")
If Not files.FileExists(interpreter) Then interpreter = "pythonw.exe"

shell.CurrentDirectory = project
command = QuotePath(interpreter) & " -X utf8 -m hischeck.gui"
For Each argument In WScript.Arguments
    command = command & " " & QuotePath(argument)
Next

On Error Resume Next
shell.Run command, 1, False
If Err.Number <> 0 Then
    MsgBox "Could not start HisCheck. Run launch_gui.cmd to view diagnostics." & vbCrLf & Err.Description, _
           vbExclamation, "HisCheck"
    WScript.Quit 1
End If

Function QuotePath(value)
    Dim trailing
    trailing = 0
    Do While trailing < Len(value)
        If Mid(value, Len(value) - trailing, 1) <> "\" Then Exit Do
        trailing = trailing + 1
    Loop
    ' Double a terminal backslash so Python receives it before the closing quote.
    QuotePath = """" & value & String(trailing, "\") & """"
End Function
