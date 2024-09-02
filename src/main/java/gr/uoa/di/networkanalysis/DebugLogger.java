package gr.uoa.di.networkanalysis;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class DebugLogger {

    private BufferedWriter writer;

    public DebugLogger(String filename) {
        try {
            writer = new BufferedWriter(new FileWriter(filename, true)); // Open in append mode
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void log(String message) {
        try {
            writer.write(message);
            writer.newLine(); // New line after each message
            writer.flush(); // Ensure data is written immediately
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void close() {
        try {
            if (writer != null) {
                writer.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
