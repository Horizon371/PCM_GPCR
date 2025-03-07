import sqlite3
import csv
import constants


def run_chembl_querry(querry, output_csv):

    # Connect to the sql database
    conn = sqlite3.connect(f"{constants.DATA_FOLDER}/chembl_35/chembl_35_sqlite/chembl_35.db")
    cursor = conn.cursor()

    # Run the database querry
    cursor.execute(querry)

    # Parse the querry results and write them to the output csv file
    results = cursor.fetchone()
    with open(output_csv, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow([i[0] for i in cursor.description])  # Write headers


    while True:
        results = cursor.fetchmany(30)
        if not results:
            break
        # Write the results to a CSV file
        with open(output_csv, 'a', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerows(results)  # Write rows

    # Close the database connection
    cursor.close()
    conn.close()