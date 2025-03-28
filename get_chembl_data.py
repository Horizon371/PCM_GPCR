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
    results = cursor.fetchall()

    with open(output_csv, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow([i[0] for i in cursor.description])  # Write headers
        csv_writer.writerows(results)

    print(f"Fetched {len(results)} rows from the database")

    # Close the database connection
    cursor.close()
    conn.close()