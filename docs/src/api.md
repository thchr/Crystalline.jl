# API

---

```@meta
CurrentModule = Crystalline
```

## Exported types
```@autodocs
Modules = [Crystalline]
Private = false
Order   = [:types]
```

## Exported methods
```@autodocs
Modules = [Crystalline]
Private = false
Order   = [:function]
```

## Unexported API

### Methods
```@autodocs
Modules = [Crystalline]
Private = true
Filter  = t->any(t′->basename(dirname(string(t′.file)))==="src", methods(t)) # restrict to methods in /src/ (e.g. exclude /build/)
```