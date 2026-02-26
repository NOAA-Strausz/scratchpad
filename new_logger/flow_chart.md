graph TD
    %% Startup Phase
    A([Start: python logger_engine.py]) --> B[Load settings.yaml]
    B --> C[(Create/Open SQLite DB)]
    C --> D[Initialize asyncio Event Loop]

    %% Task Spawning
    D --> E{Read Sensors from Config}
    E -- Type: Serial --> F[Spawn Serial Task]
    E -- Type: UDP --> G[Spawn UDP Task]
    E -- All spawned --> H((Main Async Loop))

    %% Listeners running concurrently
    F -.-> H
    G -.-> H

    subgraph Hardware & Network Interfaces
        I[Serial Listener<br>Auto-reconnects on error] -->|Awaits data| K[Read string]
        J[UDP Listener] -->|Awaits packet| K
    end

    H -.-> I
    H -.-> J

    %% The Dispatcher / Storage Pipeline
    K -->|Passes raw string| L[Dispatcher function<br>handle_incoming_data]

    subgraph Dual-Storage Pipeline
        L --> M[Generate UTC Timestamp]
        
        %% Raw Archive Path
        M --> N[(Append to Daily<br>Raw .dat File)]
        
        %% Database Path
        M --> O{yaml flag:<br>parse_data = true?}
        O -- Yes --> P[Run specific protocol parser<br>Return JSON]
        O -- No --> Q[Pass raw string directly]
        
        P --> R[(Insert Row into SQLite)]
        Q --> R
    end

    %% Graceful Shutdown Sequence
    S([User presses Ctrl+C]) -.-> T[Trigger KeyboardInterrupt]
    T --> U[Send Cancel Signal to all Tasks]
    U --> V[Wait for tasks to stop]
    V --> W[(Commit & Close SQLite DB)]
    W --> X([Clean Exit])

    classDef database fill:#f9d0c4,stroke:#333,stroke-width:2px;
    class C,N,R,W database;
